import sys
import pandas as pd
from helpers.geometry_helpers import average_molecule
from helpers.plot_functions import plot_fragment_colored

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_vdw_surface.py <path/to/inputfile>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    fragments_df = pd.read_csv(inputfilename)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    avg_fragment = average_molecule(fragments_df)
    avg_fragment.set_vdw_radii("data/vdw_radii.csv")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)
    
    # Create a sphere around the nitrogen atom
    ax, _ = plot_vdw_spheres(avg_fragment, ax)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)

    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")

    plt.show()

def plot_vdw_spheres(avg_fragment, ax):
    for atom in avg_fragment.atoms.values():
        spheres = []

        r = atom.vdw_radius

        theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

        x = r*np.sin(phi)*np.cos(theta) + atom.x
        y = r*np.sin(phi)*np.sin(theta) + atom.y
        z = r*np.cos(phi) + atom.z

        sphere = ax.plot_surface(x, y, z, color='pink', alpha=0.2, linewidth=0)

        spheres.append(sphere)
    return ax, spheres


if __name__ == "__main__":
    main()