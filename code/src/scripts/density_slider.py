import math
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import Button, Slider
from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings
from helpers.plot_functions import plot_fragment_colored, plot_density


def main():

    if len(sys.argv) != 3:
        print("Usage: python analyze_density.py <path/to/inputfile> <atom or center to count>")
        sys.exit(1)

    settings = Settings(sys.argv[1])

    # resolution of the bins, in Angstrom
    settings.set_resolution(0.25)
    settings.set_atom_to_count(sys.argv[2])

    avg_fragment = pd.read_csv(settings.get_avg_frag_filename())

    make_plot(avg_fragment, settings)


def make_plot(avg_fragment, settings):
    df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())

    df[settings.to_count_contact] = df[settings.to_count_contact] / df[settings.to_count_contact].sum()

    maximum = df[settings.to_count_contact].max()
    print(df[df[settings.to_count_contact] == maximum])

    global lower_limit, p
    lower_limit = 0.2 * maximum

    fig = plt.figure(figsize=(8, 5))

    plt.subplots_adjust(bottom=0.30)
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)
    p, ax = plot_density(ax=ax, df=df, settings=settings, lower_lim=lower_limit)

    title = f"{settings.central_group_name}--{settings.contact_group_name} 4D density plot\n\
               Resolution: {settings.resolution}"

    ax.set_title(title)

    xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    xlim, ylim, zlim = (math.floor(xlim[0]), math.ceil(xlim[1])), (math.floor(ylim[0]), math.ceil(ylim[1])),\
                       (math.floor(zlim[0]), math.ceil(zlim[1]))

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

    cax = fig.add_axes([0.8, 0.25, 0.03, 0.6])
    fig.colorbar(p, cax=cax, pad=0.2)

    axcolor = 'lightgoldenrodyellow'
    ax_resolution = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    resolution = Slider(ax_resolution, 'Res', 0.2, 1.5, valinit=0.2, valstep=0.1)

    ax_lowerlim = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    lowerlim = Slider(ax_lowerlim, 'Lim', 0, 1, valinit=0.2, valstep=0.01)

    percentage = round(df[df[settings.to_count_contact] >= lower_limit][settings.to_count_contact].sum() * 100, 2)
    text_holder = fig.text(.25, .05, f"Showing {percentage}% of all data")

    # update everything
    def update_res(val):
        print("Changed resolution to:", round(val, 5))
        settings.set_resolution(round(val, 2))
        df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())

        global lower_limit, p

        ax.cla()
        cax.cla()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)

        title = f"{settings.central_group_name}--{settings.contact_group_name} 4D density plot\n\
                   Resolution: {settings.resolution}"
        ax.set_title(title)

        p, ax1 = plot_density(ax=ax, df=df, settings=settings, lower_lim=lower_limit)
        ax1 = plot_fragment_colored(ax1, avg_fragment)

        percentage = round(df[df[settings.to_count_contact] >= lower_limit][settings.to_count_contact].sum() * 100, 2)
        text_holder.set_text(f"Showing {percentage}% of all data")

        fig.colorbar(p, cax=cax, pad=0.2)

    def update_lim(val):
        print("Changed lower limit too:", round(val, 2))
        print("Resolution:", settings.resolution)
        global lower_limit, p
        lower_limit = val * maximum

        ax.cla()
        cax.cla()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)

        _, ax1 = plot_density(ax=ax, df=df, settings=settings, lower_lim=lower_limit)
        ax1 = plot_fragment_colored(ax1, avg_fragment)

        title = f"{settings.central_group_name}--{settings.contact_group_name} 4D density plot\n\
                   Resolution: {settings.resolution}"
        ax.set_title(title)

        percentage = round(df[df[settings.to_count_contact] >= lower_limit][settings.to_count_contact].sum() * 100, 2)

        text_holder.set_text(f"Showing {percentage}% of all data")
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
