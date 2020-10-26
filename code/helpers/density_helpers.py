import math

import numpy as np
import pandas as pd


def find_maximum(df):
    """ Returns the largest coordinate for each of the x, y and z axis. """

    maxx = df.atom_x.max()
    maxy = df.atom_y.max()
    maxz = df.atom_z.max()

    return maxx, maxy, maxz


def find_minimum(df):
    """ Returns the smallest coordinate for each of the x, y and z axis. """

    minx = df.atom_x.min()
    miny = df.atom_y.min()
    minz = df.atom_z.min()

    return minx, miny, minz


def calculate_no_bins(minimum, maximum, resolution):
    bins_neg = math.ceil(abs(minimum + 0.5 * resolution)/resolution)
    bins_pos = math.ceil(abs(maximum - 0.5 * resolution)/resolution)

    minimum, maximum = -bins_neg * resolution, bins_pos * resolution

    # add one extra bin so we can put the origin in the middle of a bin
    no_bins = bins_neg + bins_pos + 1

    return no_bins, minimum, maximum


def prepare_df(fragments_df, settings):
    maxx, maxy, maxz = find_maximum(fragments_df)
    minx, miny, minz = find_minimum(fragments_df)

    no_bins_x, minx, maxx = calculate_no_bins(minx, maxx, settings.resolution)
    no_bins_y, miny, maxy = calculate_no_bins(miny, maxy, settings.resolution)
    no_bins_z, minz, maxz = calculate_no_bins(minz, maxz, settings.resolution)

    amount_bins = no_bins_x * no_bins_y * no_bins_z
    indices = [i for i in range(0, amount_bins)]

    bins = [np.linspace(minx, maxx, num=no_bins_x, endpoint=False),
            np.linspace(miny, maxy, num=no_bins_y, endpoint=False),
            np.linspace(minz, maxz, num=no_bins_z, endpoint=False)] 

    df = add_boundaries_per_bin(bins, indices, settings)

    return df


def add_boundaries_per_bin(bins, indices, settings):

    df = pd.DataFrame([], index=indices)

    bins_x, bins_y, bins_z = bins[0], bins[1], bins[2]

    xl, yl, zl = len(bins_x), len(bins_y), len(bins_z)

    xstart_list = np.repeat(bins_x, (yl * zl))
    ystart_list = list(np.repeat(bins_y, zl)) * xl
    zstart_list = list(bins_z) * (xl * yl)

    df["xstart"] = xstart_list
    df["ystart"] = ystart_list
    df["zstart"] = zstart_list

    df[settings.to_count_contact] = 0.0

    df = df.apply(pd.to_numeric, downcast='float', errors='coerce')

    return df
