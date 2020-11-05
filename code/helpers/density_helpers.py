import math

import numpy as np
import pandas as pd

from numba import jit
from numba import prange


def calculate_no_bins(minimum, maximum, resolution):
    """ Calculates the number of bins needed between a minimum and a maximum at a certain resolution.
        Returns the number of  bins, and the adjusted minimum and maximum to make the bins fit precisely. """

    bins_neg = math.ceil(abs(minimum + 0.5 * resolution)/resolution)
    bins_pos = math.ceil(abs(maximum - 0.5 * resolution)/resolution)

    minimum, maximum = -bins_neg * resolution, bins_pos * resolution

    # add one extra bin so we can put the origin in the middle of a bin
    no_bins = bins_neg + bins_pos + 1

    return no_bins, minimum, maximum


def prepare_df(df, settings):
    maxx, maxy, maxz = df.atom_x.max(), df.atom_y.max(), df.atom_z.max()
    minx, miny, minz = df.atom_x.min(), df.atom_y.min(), df.atom_z.min()

    no_bins_x, minx, maxx = calculate_no_bins(minx, maxx, settings.resolution)
    no_bins_y, miny, maxy = calculate_no_bins(miny, maxy, settings.resolution)
    no_bins_z, minz, maxz = calculate_no_bins(minz, maxz, settings.resolution)

    amount_bins = no_bins_x * no_bins_y * no_bins_z
    indices = [i for i in range(0, amount_bins)]

    bins = [np.linspace(minx, maxx, num=no_bins_x, endpoint=False),
            np.linspace(miny, maxy, num=no_bins_y, endpoint=False),
            np.linspace(minz, maxz, num=no_bins_z, endpoint=False)]

    df = add_boundaries_per_bin(bins, indices)
    df[settings.to_count_contact] = 0

    return df


def add_boundaries_per_bin(bins, indices):

    df = pd.DataFrame([], index=indices)

    bins_x, bins_y, bins_z = bins[0], bins[1], bins[2]

    xl, yl, zl = len(bins_x), len(bins_y), len(bins_z)

    xstart_list = np.repeat(bins_x, (yl * zl))
    ystart_list = list(np.repeat(bins_y, zl)) * xl
    zstart_list = list(bins_z) * (xl * yl)

    df["xstart"] = xstart_list
    df["ystart"] = ystart_list
    df["zstart"] = zstart_list

    df = df.apply(pd.to_numeric, downcast='float', errors='coerce')

    return df


@jit(nopython=True, parallel=True)
def calc_distances(in_vdw_volume, bin_coordinates, avg_f_p, indices, extra):

    for i in prange(len(indices)):
        idx = indices[i]

        bin_point = bin_coordinates[idx[0]]
        distance = np.sqrt(np.sum((bin_point - avg_f_p[:3])**2))

        if distance < avg_f_p[3] + extra:
            in_vdw_volume[idx] = 1

    return in_vdw_volume


def find_min_max_bounds(avg_fragment, extra):
    avg_fragment['minx'] = avg_fragment['atom_x'] - avg_fragment['vdw_radius'] - extra
    avg_fragment['miny'] = avg_fragment['atom_y'] - avg_fragment['vdw_radius'] - extra
    avg_fragment['minz'] = avg_fragment['atom_z'] - avg_fragment['vdw_radius'] - extra
    avg_fragment['maxx'] = avg_fragment['atom_x'] + avg_fragment['vdw_radius'] + extra
    avg_fragment['maxy'] = avg_fragment['atom_y'] + avg_fragment['vdw_radius'] + extra
    avg_fragment['maxz'] = avg_fragment['atom_z'] + avg_fragment['vdw_radius'] + extra

    return avg_fragment


def count_bins_in_vdw(avg_fragment, extra):
    """ Input: the coordinates of the average fragment, the resolution on which we are calculating and the radius
        of the contact group + 0.5.

        Output: dataframe with defined bins and whether they are in the volume or not, and the amount of bins that is
        in the vdw volume. """

    # 0.1 is accurate and semi-instantaneous calculation, see thesis
    resolution = 0.1

    avg_fragment = find_min_max_bounds(avg_fragment, extra)

    minx, miny, minz = avg_fragment['minx'].min(), avg_fragment['miny'].min(), avg_fragment['minz'].min()
    maxx, maxy, maxz = avg_fragment['maxx'].max(), avg_fragment['maxy'].max(), avg_fragment['maxz'].max()

    no_bins_x, minx, maxx = calculate_no_bins(minx, maxx, resolution)
    no_bins_y, miny, maxy = calculate_no_bins(miny, maxy, resolution)
    no_bins_z, minz, maxz = calculate_no_bins(minz, maxz, resolution)

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

        fragment_point = np.array([atom.atom_x, atom.atom_y, atom.atom_z, atom.vdw_radius])

        in_vdw_vol = calc_distances(in_vdw_vol, bin_coordinates, fragment_point, indices, extra)

    total = np.sum(in_vdw_vol)

    return total * resolution**3


def find_available_volume(avg_fragment, extra):
    volume_max = count_bins_in_vdw(avg_fragment=avg_fragment, extra=extra)
    volume_central = count_bins_in_vdw(avg_fragment=avg_fragment, extra=0)

    return volume_max - volume_central
