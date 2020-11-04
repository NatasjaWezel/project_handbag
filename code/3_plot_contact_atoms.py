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
from helpers.geometry_helpers import (get_vdw_distance_contact,
                                      make_avg_fragment_if_not_exists,
                                      make_coordinate_df)
from helpers.headers import AXCOLOR
from helpers.helpers import read_results_alignment
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres


def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_all_contact_atoms.py <path/to/inputfile> <atom/center to count>")
        sys.exit(1)

    settings = Settings(sys.argv[1])
    settings.set_atom_to_count(sys.argv[2])

    df = read_results_alignment(settings.get_aligned_csv_filename())

    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    # grab only the atoms that are in the contact groups
    df = df[df.in_central_group == False]

    vdw_distance_contact = get_vdw_distance_contact(df, settings)

    coordinate_df = make_coordinate_df(df, settings, avg_fragment)
    make_plot(avg_fragment, coordinate_df, vdw_distance_contact)


def make_plot(avg_fragment, coordinate_df, longest_vdw_contact):
    """ Plot all the surrounding contact groups around the central group. """

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    # plot the (average of the) central group
    ax = plot_fragment_colored(ax, avg_fragment)
    ax, _ = plot_vdw_spheres(avg_fragment, ax, color='pink')

    points = ax.scatter(list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom
                             + longest_vdw_contact].atom_x),
                        list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom
                             + longest_vdw_contact].atom_y),
                        list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom
                             + longest_vdw_contact].atom_z), s=1, c="red")

    vdw_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=AXCOLOR)
    vdw_slider = Slider(vdw_slider_ax, 'VDW radius + ', -2, 3, valinit=0, valstep=0.1)

    def update(val):
        val = round(val, 2)

        show_df = coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact + val]

        print(len(coordinate_df), len(show_df), val)
        points._offsets3d = (list(show_df.atom_x), list(show_df.atom_y), list(show_df.atom_z))
        fig.canvas.draw_idle()

    vdw_slider.on_changed(update)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


if __name__ == "__main__":
    main()
