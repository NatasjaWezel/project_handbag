import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    # more bins = higher resolution
    amount_bins = 10

    points = [(-4, -4), (2.2, 4.6), (2.2, 4.6), (5.4, 6.8), (2.4, 3.7), (3.0, 3.0)]

    maxx = 0
    maxy = 0

    for point in points:
        if abs(point[0]) > maxx:
            maxx = abs(point[0])
        
        if abs(point[1]) > maxy:
            maxy = abs(point[1])

    maximum = max([maxx, maxy])
    print(maximum)
    minimum = -maximum

    bins = np.linspace(minimum, maximum, num=amount_bins+1)
    print("bins: ", bins)

    indices = [i for i in range(0, amount_bins * amount_bins)]
    df = pd.DataFrame(columns=['xstart', 'xend', 'ystart', 'yend', 'amount'], index=indices)

    print(df.head())

    current_index = 0
    for i in range(0, amount_bins):
        xstart, xend = bins[i], bins[i + 1]

        for j in range(0, amount_bins):
            
            ystart, yend = bins[j], bins[j + 1]
            amount = 0

            df.loc[current_index] = [xstart, xend, ystart, yend, amount]
            current_index += 1

    df = count_points_per_square(df, points)

    plot_density(df)
    

def plot_density(df):
    df['ymiddle'] = (df['ystart'] + df['yend']) / 2
    df['xmiddle'] = (df['xstart'] + df['xend']) / 2


    # TODO: think of a 'normal' coloring way :)
    df['color'] = None
    df.loc[df.amount == 1, 'color'] = 'lightblue'
    df.loc[df.amount == 2, 'color'] = 'blue'
    
    points = df[df.amount > 0]
    
    for _, point in points.iterrows():
        plt.scatter(point.xmiddle, point.ymiddle, s=point.amount*1000, c=point.color)

    plt.show()

def count_points_per_square(df, points):
    for point in points:
        x = point[0]
        y = point[1]

        # TODO: find proper way to deal with edge cases
        print(x,y, len(df[(df.xstart < x) & (df.xend >= x) & (df.ystart < y) & (df.yend >= y)]))
        
        df.loc[(df.xstart < x) & (df.xend >= x) & (df.ystart < y) & (df.yend >= y), 'amount'] = df[(df.xstart < x) & (df.xend >= x) & (df.ystart < y) & (df.yend >= y)].amount + 1
    
    return df

def plot_original_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for point in points:
        ax.scatter(point[0], point[1])

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
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