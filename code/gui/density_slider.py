import sys

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import math

from matplotlib.widgets import Slider, Button

from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings
from calc_avg_fragment_2 import make_avg_fragment_if_not_exists, read_results_alignment
from helpers.plot_functions import plot_density, plot_fragment_colored


def main():

    if len(sys.argv) != 3:
        print("Usage: python analyze_density.py <path/to/inputfile> <atom or center to count>")
        sys.exit(1)

    settings = Settings(sys.argv[1])

    # resolution of the bins, in Angstrom
    settings.set_resolution(0.5)
    settings.set_atom_to_count(sys.argv[2])

    df = read_results_alignment(settings.get_aligned_csv_filename())
    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    make_plot(avg_fragment, settings)

def get_density_data(settings):
    # first density: 0.2
    df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())

    df['ymiddle'] = (df['ystart'] * 2 + settings.resolution) / 2
    df['xmiddle'] = (df['xstart'] * 2 + settings.resolution) / 2
    df['zmiddle'] = (df['zstart'] * 2 + settings.resolution) / 2

    # normalize per column
    df[settings.to_count_contact] = df[settings.to_count_contact] / df[settings.to_count_contact].sum()
    
    return df


def make_plot(avg_fragment, settings):

    lower_limit = 0.00001
    df = get_density_data(settings)
    points = df[df[settings.to_count_contact] > lower_limit]
    min_low, max_low = df[settings.to_count_contact].min(), df[settings.to_count_contact].max()

    fig = plt.figure(figsize=(8,5))
    plt.subplots_adjust(bottom=0.30)
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    norm = plt.Normalize(lower_limit, points[settings.to_count_contact].max())
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue", "fuchsia", "red"])

    p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle),
                   s=list(10000 * points[settings.to_count_contact]),
                   c=list(points[settings.to_count_contact]),
                   cmap=cmap,
                   norm=norm)

    ax.set_title("4D density plot\n Resolution: " + str(settings.resolution))
    xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    xlim, ylim, zlim = (math.floor(xlim[0]), math.ceil(xlim[1])), (math.floor(ylim[0]), math.ceil(ylim[1])), (math.floor(zlim[0]), math.ceil(zlim[1]))

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

    cax = fig.add_axes([0.8, 0.25, 0.03, 0.6])
    colorbar = fig.colorbar(p, cax=cax, pad=0.2)

    axcolor = 'lightgoldenrodyellow'
    ax_resolution = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    resolution = Slider(ax_resolution, 'Res', 0.2, 1.5, valinit=0.5, valstep=0.1)

    ax_lowerlim = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    lowerlim = Slider(ax_lowerlim, 'Lim', min_low, max_low, valinit=0.00001, valstep=0.00001)

    # update everything
    def update_res(val):
        print("Changed resolution to:", round(val, 5))
        settings.set_resolution(round(val, 2))
        df = get_density_data(settings)
        points = df[df[settings.to_count_contact] > lower_limit]
        ax.cla()
        cax.cla()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)

        norm = plt.Normalize(lower_limit, points[settings.to_count_contact].max())
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue", "fuchsia", "red"])

        p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle),
                    s=list(10000 * points[settings.to_count_contact]),
                    c=list(points[settings.to_count_contact]),
                    cmap=cmap,
                    norm=norm)

        plot_fragment_colored(ax, avg_fragment)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        ax.set_title("4D density plot\n Resolution: " + str(settings.resolution))

        fig.colorbar(p, cax=cax, pad=0.2)

    def update_lim(val):
        print("Changed lower limit too:", round(val, 2))
        lower_limit = val
        points = df[df[settings.to_count_contact] > val]

        ax.cla()
        cax.cla()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)

        norm = plt.Normalize(lower_limit, points[settings.to_count_contact].max())
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue", "fuchsia", "red"])

        p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle),
                    s=list(10000 * points[settings.to_count_contact]),
                    c=list(points[settings.to_count_contact]),
                    cmap=cmap,
                    norm=norm)

        plot_fragment_colored(ax, avg_fragment)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        ax.set_title("4D density plot\n Resolution: " + str(settings.resolution))

        fig.colorbar(p, cax=cax, pad=0.2)


    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


    def reset(event):
        resolution.reset()
        lowerlim.reset()

    button.on_clicked(reset)
    resolution.on_changed(update_res)
    lowerlim.on_changed(update_lim) 

    plt.show()


if __name__ == "__main__":
    main()
