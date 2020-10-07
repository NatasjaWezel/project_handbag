import sys
import pandas as pd
from helpers.geometry_helpers import make_avg_fragment_if_not_exists
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres
from helpers.helpers import read_results_alignment

from classes.Settings import Settings

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.widgets import RadioButtons, CheckButtons, Slider

import numpy as np

from helpers.headers import AXCOLOR

def main():
    if len(sys.argv) != 3:
        print("Usage: python plot_vdw_surface.py <path/to/inputfile>")
        sys.exit(1)

    settings = Settings(sys.argv[1])
    settings.set_central_group()

    df = read_results_alignment(settings.get_aligned_csv_filename())

    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    # to plot the vdw surface we need the vdw radii
    for atom in avg_fragment.atoms.values():
        radius = settings.get_vdw_radius(atom.symbol)
        print(atom, radius)
        atom.set_vdw_radius(radius)

    plot_vdw_surface(avg_fragment)

def plot_vdw_surface(avg_fragment):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.25, bottom=0.25)

    ax.margins(x=0)
    ax = plot_fragment_colored(ax, avg_fragment)
    # Create vdw spheres around central group atoms
    ax, spheres = plot_vdw_spheres(avg_fragment, ax)

    # visibility of the spheres    
    check_ax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=AXCOLOR)
    radio = CheckButtons(check_ax, ['visible'])

    def switch_visibility(label):
        for sphere in spheres:
            sphere.set_visible(not sphere.get_visible())

        fig.canvas.draw_idle()

    radio.on_clicked(switch_visibility)
    
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")

    plt.show()


if __name__ == "__main__":
    main()