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
    
    # Create a sphere
    r = 1

    phi, theta = np.mgrid[0.0: np.pi : 100j, 0.0 : 2.0*np.pi : 100j]
    x = r*np.sin(phi)*np.cos(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(phi)

    ax.plot_surface(x, y, z, color='red', alpha=0.3, linewidth=0)

    plt.show()


if __name__ == "__main__":
    main()