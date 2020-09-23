import sys
import pandas as pd
from helpers.geometry_helpers import average_fragment
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres
from helpers.helpers import read_results_alignment

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.widgets import RadioButtons, CheckButtons, Slider

import numpy as np

from helpers.headers import AXCOLOR

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_vdw_surface.py <path/to/inputfile>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    prefix = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] 
    avg_fragment_name = prefix + "_avg_fragment.pkl"

    df = read_results_alignment(inputfilename)

    avg_fragment = average_fragment(avg_fragment_name, df)

    # to plot the vdw surface we need the vdw radii
    avg_fragment.set_vdw_radii("data/vdw_radii.csv")

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