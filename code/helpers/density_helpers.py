import numpy as np
import pandas as pd
import math


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
   

def prepare_df(fragments_df, resolution):
    maxx, maxy, maxz = find_maximum(fragments_df)
    minx, miny, minz = find_minimum(fragments_df)

    # TODO: check if this is entirely correct
    # i think the volumes still differ per search now
    no_bins_x = math.ceil((maxx + abs(minx))/resolution)
    no_bins_y = math.ceil((maxy + abs(miny))/resolution)
    no_bins_z = math.ceil((maxz + abs(minz))/resolution)

    amount_bins = no_bins_x * no_bins_y * no_bins_z
    indices = [i for i in range(0, amount_bins)]

    bins = [np.linspace(minx, maxx, num=no_bins_x + 1), np.linspace(miny, maxy, num=no_bins_y + 1), np.linspace(minz, maxz, num=no_bins_z + 1)] 
    
    volume_total = abs(maxx - minx) * abs(maxy - miny) * abs(maxz - minz)
    print("Volume per bin & volume per bin according to resolution: ", resolution**3, volume_total/amount_bins)

    df = add_boundaries_per_bin(bins, indices)
   
    return df


def add_boundaries_per_bin(bins, indices):
    
    df = pd.DataFrame([], index=indices)

    bins_x, bins_y, bins_z = bins[0], bins[1], bins[2]
    
    # TODO: fix -1 earlier
    xl = len(bins_x) -1
    yl = len(bins_y) -1 
    zl = len(bins_z) -1
    
    print(xl, yl, zl)

    xstart_list = np.repeat(bins_x[:-1], (yl * zl))
    ystart_list = list(np.repeat(bins_y[:-1], zl)) * xl 
    zstart_list = list(bins_z[:-1]) * (xl * yl)

    df["xstart"] = xstart_list
    df["ystart"] = ystart_list
    df["zstart"] = zstart_list

    df = df.apply(pd.to_numeric, downcast='float', errors='coerce')

    return df


def add_one_to_bin(df, column_name, coordinates, resolution):
    x, y, z = coordinates[0], coordinates[1], coordinates[2]

    # TODO: find out what happens with points exactly on a bin-line
    index = df.index[((df.xstart <= x) & (df.xstart + resolution >= x) & 
                        (df.ystart <= y) & (df.ystart + resolution >= y) & 
                        (df.zstart <= z) & (df.zstart + resolution >= z))]

    df.loc[index, column_name] = df.loc[index, column_name] + 1

    return df

