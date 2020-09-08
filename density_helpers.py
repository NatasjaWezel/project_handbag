import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():
    # more bins = higher resolution
    amount_bins = 10

    points = [(-4, -4, 2), (2.2, 4.6, 3), (2.2, 4.6, 5.3), (5.4, 6.8, 5.1), (5.3, 6.8, 5.1), (2.4, 3.7, -2.6),
                        (5.3, 6.8, 5.1), (2.4, 3.7, -2.6), (-4, -6, 2.3),
                        (3.0, 3.0, 3.0), (3.0, 3.0, 3.0), (3.0, 3.2, 3.0),
                        (3.0, 3.0, 3.0), (3.0, 3.0, 3.0), (3.0, 3.2, 3.0),
                        (-7.9, -5.6, -6.5)]

    maximum = find_maximum(points)
    minimum = -maximum
    print(minimum, maximum)

    bins = np.linspace(minimum, maximum, num=amount_bins+1)
    indices = [i for i in range(0, amount_bins * amount_bins * amount_bins)]

    empty_df = prepare_df(amount_bins, bins, indices)

    df = count_points_per_square(empty_df, points)
    test_count(points, df)
    
    plot_density(df, amount_bins, minimum, maximum)

def test_count(points, df):
    assert (len(points) == df.amount.sum()), "Amount of points is not divided into bins correctly"

def find_maximum(points):
    maxx, maxy, maxz = 0, 0, 0

    for point in points:
        if abs(point[0]) > maxx:
            maxx = abs(point[0])
        
        if abs(point[1]) > maxy:
            maxy = abs(point[1])

        if abs(point[2] > maxz):
            maxz = abs(point[2])

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
        p = ax.scatter(point.xmiddle, point.ymiddle, point.zmiddle, s=point.amount*250, c=point.amount, cmap=cmap, norm=norm)

    ax.set_title("4D density plot\n Bins: " + str(amount_bins))
    ax.set_xlim(minimum, maximum)
    ax.set_ylim(minimum, maximum)
    ax.set_zlim(minimum, maximum)

    fig.colorbar(p)
    plt.show()

def count_points_per_square(df, points):
    for point in points:
        x, y, z = point[0], point[1], point[2]
        
        # TODO: find out what happens with points exactly on a bin-line
        df.loc[(df.xstart <= x) & (df.xend >= x) & (df.ystart <= y) & (df.yend >= y) & (df.zstart <= z) & (df.zend >= z),
                    'amount'] = df[(df.xstart <= x) & (df.xend >= x) & (df.ystart <= y) & (df.yend >= y) & (df.zstart <= z) & (df.zend >= z)].amount + 1

    return df

def plot_original_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for point in points:
        ax.scatter(point[0], point[1])

    # Move left y-axis and bottom x-axis to centre, passing through (0,0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)

    plt.savefig("points.png")


if __name__ == "__main__":
    main()