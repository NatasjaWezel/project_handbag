# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments, and then plots the central
# group and all contact atoms/the centers of the contact groups around it.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings
from classes.Radii import Radii
from helpers.geometry_helpers import (make_coordinate_df)
from constants.colors import AXCOLOR
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres

from constants.paths import RADII_CSV, WORKDIR

import pandas as pd


def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_all_contact_atoms.py <path/to/inputfile> <atom/center to count>")
        sys.exit(1)

    settings = Settings(WORKDIR, sys.argv[1])
    settings.set_atom_to_count(sys.argv[2])

    df = pd.read_csv(settings.get_aligned_csv_filename())

    avg_fragment = pd.read_csv(settings.get_avg_frag_filename())

    # grab only the atoms that are in the contact groups
    df = df[df.label == "-"]

    radii = Radii(RADII_CSV)
    vdw_distance_contact = radii.get_vdw_distance_contact(settings.to_count_contact)

    coordinate_df = make_coordinate_df(df, settings, avg_fragment, radii)
    print(coordinate_df.x.max(), coordinate_df.y.max(), coordinate_df.z.max())
    print(coordinate_df[coordinate_df.x == coordinate_df.x.max()])

    title = "Central fragment: " + settings.central_name + "\n" +\
            "Scattered contact atoms: " + settings.contact_name + "(" + settings.to_count_contact + ")"
    make_plot(avg_fragment, coordinate_df, vdw_distance_contact, title)


def make_plot(avg_fragment, coordinate_df, longest_vdw_contact, title):
    """ Plot all the surrounding contact groups around the central group. """

    # calculate corrected vdw distance
    coordinate_df["vdw_corr"] = coordinate_df.distance - coordinate_df.vdw_closest_atom - longest_vdw_contact

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(bottom=0.25)
    plt.title(title)

    # plot the (average of the) central group
    ax = plot_fragment_colored(ax, avg_fragment)

    # TODO: make checkbox
    # ax, _ = plot_vdw_spheres(avg_fragment, ax)

    xlim, ylim, zlim = list(ax.get_xlim()), list(ax.get_ylim()), list(ax.get_zlim())
    minn = min([xlim[0], ylim[0], zlim[0]])
    maxx = max([xlim[1], ylim[1], zlim[1]])

    ax.set_xlim((minn, maxx))
    ax.set_ylim((minn, maxx))
    ax.set_zlim((minn, maxx))

    print(coordinate_df.vdw_corr.min(), coordinate_df.vdw_corr.max())

    points = ax.scatter(list(coordinate_df[coordinate_df.vdw_corr <= 0.5].x),
                        list(coordinate_df[coordinate_df.vdw_corr <= 0.5].y),
                        list(coordinate_df[coordinate_df.vdw_corr <= 0.5].z), s=1, c="red")

    vdw_slider_ax = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=AXCOLOR)
    vdw_slider = Slider(vdw_slider_ax, 'Upper limit', round(coordinate_df.vdw_corr.min(), 2), 0.5,
                        valinit=0.5, valstep=0.01)

    vdw_slider_ax_lo = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=AXCOLOR)
    vdw_slider_lo = Slider(vdw_slider_ax_lo, 'Lower limit', round(coordinate_df.vdw_corr.min(), 2), 0.5,
                           valinit=round(coordinate_df.vdw_corr.min(), 2), valstep=0.01)

    vdw_slider_ax_lo.text(x=-2.15, y=1.5, s="Show contacts in van der Waals (vdW) corrected distance")

    def update(val):
        lo = vdw_slider_lo.val
        hi = vdw_slider.val

        show_df = coordinate_df[(coordinate_df.vdw_corr >= lo) &
                                (coordinate_df.vdw_corr <= hi)]

        print(f"{hi :2f}, {lo :2f}, {len(coordinate_df)}, {len(show_df)}")

        points._offsets3d = (list(show_df.x), list(show_df.y), list(show_df.z))
        fig.canvas.draw_idle()

    vdw_slider.on_changed(update)
    vdw_slider_lo.on_changed(update)

    ax.set_xlabel('X axis ($\\AA$)')
    ax.set_ylabel('Y axis ($\\AA$)')
    ax.set_zlabel('Z axis ($\\AA$)')

    plt.show()


if __name__ == "__main__":
    main()
