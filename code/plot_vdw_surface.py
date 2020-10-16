import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings
from helpers.geometry_helpers import (calculate_longest_vdw_radius_contact,
                                      make_avg_fragment_if_not_exists)
from helpers.headers import AXCOLOR
from helpers.helpers import read_results_alignment
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_vdw_surface.py <path/to/inputfile>")
        sys.exit(1)

    settings = Settings(sys.argv[1])

    df = read_results_alignment(settings.get_aligned_csv_filename())

    first_contact_group = df[(df.id == df.id.unique()[0]) & (df.in_central_group == False)]
    longest_vdw = calculate_longest_vdw_radius_contact(first_contact_group, settings)

    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    # to plot the vdw surface we need the vdw radii
    for atom1 in avg_fragment.atoms.values():
        radius = settings.get_vdw_radius(atom1.symbol)
        atom1.set_vdw_radius(radius)

    plot_vdw_surface(avg_fragment, longest_vdw)


def plot_vdw_surface(avg_fragment, longest_vdw_contact):
    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.25, bottom=0.25)

    ax.margins(x=0)
    ax = plot_fragment_colored(ax, avg_fragment)

    # Create vdw spheres around central group atoms
    ax, spheres1 = plot_vdw_spheres(avg_fragment, ax, 'red')

    # Create vdw spheres around central group atoms
    ax, spheres2 = plot_vdw_spheres(avg_fragment, ax, 'orange', extra=0.5)

    # Create vdw spheres around central group atoms
    ax, spheres3 = plot_vdw_spheres(avg_fragment, ax, 'pink', extra=0.5+longest_vdw_contact)

    # visibility of the spheres
    check_ax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=AXCOLOR)
    radio = CheckButtons(check_ax, ['vdwradius', 'vdwradis + 0.5', 'vdwradius + 0.5 + contact vdw'])

    def switch_visibility(label):
        if label == 'vdwradius':
            for sphere in spheres1:
                sphere.set_visible(not sphere.get_visible())
        elif label == 'vdwradis + 0.5':
            for sphere in spheres2:
                sphere.set_visible(not sphere.get_visible())
        else:
            for sphere in spheres3:
                sphere.set_visible(not sphere.get_visible())

        fig.canvas.draw_idle()

    radio.on_clicked(switch_visibility)

    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")

    plt.show()


if __name__ == "__main__":
    main()
