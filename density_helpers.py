import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def test_count(density_df, points_df):
    assert (len(points_df) == density_df.amount.sum()), "Amount of points is not divided into bins correctly"

def find_maximum(df):

    maxx = df.atom_x.abs().max()
    maxy = df.atom_y.abs().max()
    maxz = df.atom_z.abs().max()

    maximum = max([maxx, maxy, maxz])
    
    return maximum

def prepare_df(amount_bins, bins, indices):
    df = pd.DataFrame(columns=['xstart', 'xend', 'ystart', 'yend', 'zstart', 'zend', 'amount'], index=indices)

    current_index = 0
    for i in range(0, amount_bins):
        xstart, xend = bins[i], bins[i + 1]

        for j in range(0, amount_bins):
            ystart, yend = bins[j], bins[j + 1]

            for k in range(0, amount_bins):
                zstart, zend = bins[k], bins[k + 1]
                
                amount = 0

                df.loc[current_index] = [xstart, xend, ystart, yend, zstart, zend, amount]
                current_index += 1
    
    return df

def plot_density(df, amount_bins, minimum, maximum):
    df['ymiddle'] = (df['ystart'] + df['yend']) / 2
    df['xmiddle'] = (df['xstart'] + df['xend']) / 2
    df['zmiddle'] = (df['zstart'] + df['zend']) / 2

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    norm = plt.Normalize(1, df.amount.max())
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue","violet","red"])
        
    points = df[df.amount > 0]
    
    for _, point in points.iterrows():
        p = ax.scatter(point.xmiddle, point.ymiddle, point.zmiddle, s=point.amount, c=point.amount, cmap=cmap, norm=norm)

    ax.set_title("4D density plot\n Bins: " + str(amount_bins))
    ax.set_xlim(minimum, maximum)
    ax.set_ylim(minimum, maximum)
    ax.set_zlim(minimum, maximum)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    fig.colorbar(p)
    plt.show()

def count_points_per_square(df, points_df):
    total = len(points_df)
    for index, point in points_df.iterrows():
        x, y, z = point.atom_x, point.atom_y, point.atom_z
        
        # TODO: find out what happens with points exactly on a bin-line
        df.loc[(df.xstart <= x) & (df.xend >= x) & (df.ystart <= y) & (df.yend >= y) & (df.zstart <= z) & (df.zend >= z),
                    'amount'] = df[(df.xstart <= x) & (df.xend >= x) & (df.ystart <= y) & (df.yend >= y) & (df.zstart <= z) & (df.zend >= z)].amount + 1

        if index % 100 == 0:
            print(index, "/", total)
        
        if index == 1000:
            break

    # test_count(df, points_df)
    return df