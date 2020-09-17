import numpy as np
import pandas as pd
import math

from helpers.tests import test_count

from helpers.Atom import Atom
from helpers.Fragment import Fragment

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
   
def prepare_df(fragments_df, resolution, to_count):
    maxx, maxy, maxz = find_maximum(fragments_df)
    minx, miny, minz = find_minimum(fragments_df)

    # TODO: check if this is entirely correct
    # i think the volumes still differ per search now
    no_bins_x = math.ceil((maxx + abs(minx))/resolution)
    no_bins_y = math.ceil((maxy + abs(miny))/resolution)
    no_bins_z = math.ceil((maxz + abs(minz))/resolution)

    amount_bins = no_bins_x * no_bins_y * no_bins_z
    indices = [i for i in range(0, amount_bins)]

    print("Bins: ", no_bins_x, no_bins_y, no_bins_z)

    bins = [np.linspace(minx, maxx, num=no_bins_x + 1),
            np.linspace(miny, maxy, num=no_bins_y + 1),
            np.linspace(minz, maxz, num=no_bins_z + 1)] 
 
    df = add_boundaries_per_bin(bins, indices)

    for column in to_count:
        df["amount_" + column] = 0
    
    return df

def add_boundaries_per_bin(bins, indices):
    
    df = pd.DataFrame(columns=['xstart', 'xend', 'ystart', 'yend', 'zstart', 'zend'], index=indices)

    current_index = 0

    bins_x, bins_y, bins_z = bins[0], bins[1], bins[2]

    for i in range(0, len(bins_x) - 1):
        xstart, xend = bins_x[i], bins_x[i + 1]

        for j in range(0, len(bins_y) - 1):
            ystart, yend = bins_y[j], bins_y[j + 1]

            for k in range(0, len(bins_z) - 1):
                zstart, zend = bins_z[k], bins_z[k + 1]
                
                df.loc[current_index] = [xstart, xend, ystart, yend, zstart, zend]
                current_index += 1
    
    return df


def add_one_to_bin(df, columname, x, y, z):
    # TODO: find out what happens with points exactly on a bin-line
    df.loc[(df.xstart <= x) & (df.xend >= x) & 
                (df.ystart <= y) & (df.yend >= y) & 
                (df.zstart <= z) & (df.zend >= z),
                columname] = df[(df.xstart <= x) & (df.xend >= x) & 
                                (df.ystart <= y) & (df.yend >= y) & 
                                (df.zstart <= z) & (df.zend >= z)][columname] + 1

    return df

