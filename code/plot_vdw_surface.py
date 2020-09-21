import sys
import pandas as pd
from helpers.geometry_helpers import average_fragment
from helpers.plot_functions import plot_fragment_colored

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

    fragments_df = pd.read_csv(inputfilename)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "atom_symbol", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    avg_fragment = average_fragment(avg_fragment_name, fragments_df)

    # to plot the vdw surface we need the vdw radii
    avg_fragment.set_vdw_radii("data/vdw_radii.csv")

    plot_vdw_surface(avg_fragment)

def plot_vdw_surface(avg_fragment):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.25, bottom=0.25)

    ax.margins(x=0)
    ax = plot_fragment_colored(ax, avg_fragment)
    # Create a sphere around the nitrogen atom
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


def plot_vdw_spheres(avg_fragment, ax, slider=0.0):
    spheres = []

    for atom in avg_fragment.atoms.values():
        r = atom.vdw_radius + slider

        theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

        x = r * np.sin(phi) * np.cos(theta) + atom.x
        y = r * np.sin(phi) * np.sin(theta) + atom.y
        z = r * np.cos(phi) + atom.z

        sphere = ax.plot_surface(x, y, z, color='pink', alpha=0.2, linewidth=0)
        spheres.append(sphere)

    assert len(spheres) == len(avg_fragment.atoms.keys()), "Something went wrong with plotting the vdw surfaces"

    return ax, spheres


if __name__ == "__main__":
    main()