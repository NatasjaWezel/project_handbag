# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `density_helpers` contains functions for making bins and calculating the amount of datapoints in the bins
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import math

import numpy as np
import pandas as pd

from numba import jit
from numba import prange


def calculate_no_bins(resolution, limits):
    """ Calculates the number of bins needed between a minimum and a maximum at a certain resolution.
        Returns the number of  bins, and the adjusted minimum and maximum to make the bins fit precisely. """

    minimum, maximum = limits[0], limits[1]
    bins_neg = math.ceil(abs(minimum)/resolution)
    bins_pos = math.ceil(abs(maximum)/resolution)

    minimum, maximum = -bins_neg * resolution - 0.5 * resolution, bins_pos * resolution + 0.5 * resolution

    # add one extra bin so we can put the origin in the middle of a bin
    no_bins = bins_neg + bins_pos + 1

    return no_bins, minimum, maximum


def make_density_df(settings, coordinate_df, again=False):
    """ Make density df if it doesn't already exists. """

    try:
        if again:
            raise KeyError
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
        print("Density df already existed, loaded from file")
    except (FileNotFoundError, KeyError):

        # prepare empty df
        empty_density_df = prepare_df(df=coordinate_df, settings=settings)

        # count data per bin
        density_df = count_data_points_per_bin(df=empty_density_df, contact_points_df=coordinate_df, settings=settings)

        # normalize
        density_df['datafrac_normalized'] = (density_df[settings.contact_rp] /
                                             density_df[settings.contact_rp].sum())

        # save so we can use the data but only change the plot - saves time :)
        density_df.to_hdf(settings.get_density_df_filename(), settings.get_density_df_key())

    return density_df


def count_data_points_per_bin(df, contact_points_df, settings):
    """ Counts the amount of datapoints that is in a bin. """

    contact_points_df = contact_points_df

    print("Counting points per bin: ")
    # prepare vector that will contain the amount
    amount = np.zeros(len(df))

    bin_coordinates = np.array([df.xstart, df.ystart, df.zstart])
    contact_coordinates = np.transpose(np.array([contact_points_df.x,
                                                 contact_points_df.y,
                                                 contact_points_df.z]))

    amount = fill_bins(amount, bin_coordinates, contact_coordinates, settings.resolution)

    assert amount.sum() == len(contact_points_df), "Something went wrong with filling bins" + str(amount.sum())\
        + " " + str(len(contact_points_df))

    df[settings.contact_rp] = amount

    return df


def fill_bins(amount, bin_coordinates, contact_coordinates, resolution):
    """ Count how many datapoints there are in each bin. """

    x, y, z = 0, 1, 2
    i = 0
    total = len(contact_coordinates)

    for i in range(total):
        cor = contact_coordinates[i]
        cor = np.array([cor[x], cor[y], cor[z]], dtype='float32')
        cor_res = np.array([cor[x] - resolution,
                            cor[y] - resolution,
                            cor[z] - resolution], dtype='float32')

        idx = np.where((bin_coordinates[x] <= cor[x]) & (cor_res[x] <= bin_coordinates[x]) &
                       (bin_coordinates[y] <= cor[y]) & (cor_res[y] <= bin_coordinates[y]) &
                       (bin_coordinates[z] <= cor[z]) & (cor_res[z] <= bin_coordinates[z]))

        if len(idx[0]) > 1:
            #  because of floating point imprecision, sometimes bins overlap, so get first bin
            amount[idx[0][0]] += 1
        elif len(idx[0]) == 0:
            print("Oh no this is going wrong")
            print(cor, cor_res)
        else:
            amount[idx] += 1

        if i % 5000 == 0:
            print(i, "/", total)

    return amount


def prepare_df(df, settings):
    """ Prepares an empty df containing the coordinates for all the bins. """

    maxx, maxy, maxz = df.x.max(), df.y.max(), df.z.max()
    minx, miny, minz = df.x.min(), df.y.min(), df.z.min()

    no_bins_x, minx, maxx = calculate_no_bins(resolution=settings.resolution, limits=[minx, maxx])
    no_bins_y, miny, maxy = calculate_no_bins(resolution=settings.resolution, limits=[miny, maxy])
    no_bins_z, minz, maxz = calculate_no_bins(resolution=settings.resolution, limits=[minz, maxz])

    amount_bins = no_bins_x * no_bins_y * no_bins_z

    print(f"Preparing df, amount of bins: {amount_bins}")
    print(f"X Range: ({minx :.2f},{maxx :.2f}), Y Range: ({miny :.2f},{maxy :.2f}), Z Range: ({minz :.2f},{maxz :.2f})")

    indices = [i for i in range(0, amount_bins)]

    bins = [np.linspace(minx, maxx, num=no_bins_x, endpoint=False),
            np.linspace(miny, maxy, num=no_bins_y, endpoint=False),
            np.linspace(minz, maxz, num=no_bins_z, endpoint=False)]

    df = add_boundaries_per_bin(bins, indices)
    df[settings.contact_rp] = 0

    return df


def add_boundaries_per_bin(bins, indices):
    """ Adds bins with boundaries to a dataframe. """

    bins_x, bins_y, bins_z = bins[0], bins[1], bins[2]

    xl, yl, zl = len(bins_x), len(bins_y), len(bins_z)

    xstart_list = np.repeat(bins_x, (yl * zl))
    ystart_list = list(np.repeat(bins_y, zl)) * xl
    zstart_list = list(bins_z) * (xl * yl)

    df = pd.DataFrame(np.array([xstart_list, ystart_list, zstart_list]).T,
                      columns=['xstart', 'ystart', 'zstart'],
                      index=indices)

    df = df.apply(pd.to_numeric, downcast='float', errors='coerce')

    return df


@jit(nopython=True, parallel=True)
def calc_distances(in_vdw_volume, bin_coordinates, avg_f_p, indices, extra):
    """ Calc distances from contact rp's to closest atom from the central group model. """

    for i in prange(len(indices)):
        idx = indices[i]

        bin_point = bin_coordinates[idx[0]]
        distance = np.sqrt(np.sum((bin_point - avg_f_p[:3])**2))

        if distance < avg_f_p[3] + extra:
            in_vdw_volume[idx] = 1

    return in_vdw_volume


def find_min_max_bounds(avg_fragment, extra):
    """ Find minimum and maximum coordinates of the central model. """

    avg_fragment['minx'] = avg_fragment['x'] - avg_fragment['vdw_radius'] - extra
    avg_fragment['miny'] = avg_fragment['y'] - avg_fragment['vdw_radius'] - extra
    avg_fragment['minz'] = avg_fragment['z'] - avg_fragment['vdw_radius'] - extra
    avg_fragment['maxx'] = avg_fragment['x'] + avg_fragment['vdw_radius'] + extra
    avg_fragment['maxy'] = avg_fragment['y'] + avg_fragment['vdw_radius'] + extra
    avg_fragment['maxz'] = avg_fragment['z'] + avg_fragment['vdw_radius'] + extra

    return avg_fragment


def calc_vdw_vol_central(avg_fragment, extra, resolution):
    """ Input: the coordinates of the average fragment, the resolution on which we are calculating and the radius
        of the contact group + 0.5.

        Output: dataframe with defined bins and whether they are in the volume or not, and the amount of bins that is
        in the vdw volume.

        0.1 provides an accurate and semi-instantaneous calculation, see thesis. """

    avg_fragment = find_min_max_bounds(avg_fragment, extra)

    minx, miny, minz = avg_fragment['minx'].min(), avg_fragment['miny'].min(), avg_fragment['minz'].min()
    maxx, maxy, maxz = avg_fragment['maxx'].max(), avg_fragment['maxy'].max(), avg_fragment['maxz'].max()

    no_bins_x, minx, maxx = calculate_no_bins(resolution=resolution, limits=[minx, maxx])
    no_bins_y, miny, maxy = calculate_no_bins(resolution=resolution, limits=[miny, maxy])
    no_bins_z, minz, maxz = calculate_no_bins(resolution=resolution, limits=[minz, maxz])

    amount_bins = no_bins_x * no_bins_y * no_bins_z
    indices = [i for i in range(0, amount_bins)]

    bins = [np.linspace(minx, maxx, num=no_bins_x, endpoint=False),
            np.linspace(miny, maxy, num=no_bins_y, endpoint=False),
            np.linspace(minz, maxz, num=no_bins_z, endpoint=False)]

    df = add_boundaries_per_bin(bins, indices)

    # the bin counts as inside the vdw radius only if the center is in that radius
    df['x_center'] = df.xstart + 0.5 * resolution
    df['y_center'] = df.ystart + 0.5 * resolution
    df['z_center'] = df.zstart + 0.5 * resolution

    part_density = df[(df.x_center < maxx) & (df.x_center > minx) &
                      (df.y_center < maxy) & (df.y_center > miny) &
                      (df.z_center < maxz) & (df.z_center > minz)]

    bin_coordinates = np.transpose(np.array([part_density.x_center, part_density.y_center, part_density.z_center]))

    in_vdw_vol = np.zeros(len(df))

    for i, atom in avg_fragment.iterrows():
        indices = np.transpose(np.where(in_vdw_vol == 0))

        fragment_point = np.array([atom.x, atom.y, atom.z, atom.vdw_radius])

        in_vdw_vol = calc_distances(in_vdw_vol, bin_coordinates, fragment_point, indices, extra)

    total = np.sum(in_vdw_vol)

    return total * resolution**3


def find_available_volume(avg_fragment, extra, total=False, resolution=0.1):
    """ Find the available volume. """

    avg_fragment["label"] = avg_fragment["label"].str.upper()
    avg_fragment_without_R = avg_fragment[~avg_fragment.label.str.contains("R")].copy()
    only_R = avg_fragment[avg_fragment.label.str.contains("R")].copy()

    if total:
        return calc_vdw_vol_central(avg_fragment=avg_fragment, extra=extra, resolution=resolution)

    volume_central = calc_vdw_vol_central(avg_fragment=avg_fragment, extra=0, resolution=resolution)
    volume_max = calc_vdw_vol_central(avg_fragment=avg_fragment_without_R, extra=extra, resolution=resolution)

    volume_R_min = 0

    if len(only_R) > 0:
        volume_R_min = calc_vdw_vol_central(avg_fragment=only_R, extra=0, resolution=resolution)

    return (volume_max) - (volume_central - volume_R_min/2)
