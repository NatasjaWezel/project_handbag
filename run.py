from headers import *

import numpy as np

import copy as copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from rotation_helpers import rotate_group

from tests import check_new_group_alignment

from helpers import load_molecule

def main():

    filenames = ["data/NO3_CO_vdw5/ABOKEJ.CO_NO3_vdw5.cif", "data/NO3_CO_vdw5/AJOWIG.CO_NO3_vdw5.cif"]
    bonds = [BONDS_ABOKEJ, BONDS_AJOWIG]

    for bond, filename in zip(bonds, filenames):
        molecule = load_molecule(filename=filename)

        # center on N atom
        atom_to_center = "N1"

        # for group in groups
        group = molecule.groups[0]
        group.center_coordinates(atom_to_center=atom_to_center)
        
        # TODO: abstract these from file so you only have to say "center on N atom"
        # now it can give a key error
        atoms_to_put_in_plane = ["O1", "O2"]
        group = perform_rotations(group, atoms_to_put_in_plane, plot=False, bonds=bond)

        # test if everything went right (if the math is allright)
        print(molecule.label, end=": ")
        check_new_group_alignment(group, atom_to_center, atoms_to_put_in_plane)



def perform_rotations(group, atoms_to_put_in_plane, plot, bonds):
    """ Performs three rotations to lie three of the atoms in the xy plane, one of those
        on the x-axis. """ 

    print("Original coordinates")
    print(group)
    
    groups_to_plot = []
    labels = []
    # molecules_to_plot.append(copy.deepcopy(molecule))
    # labels.append("original")

    """ First rotation: puts first atom on xy-plane if it already was on a plane, 
        and above the x-axis if it wasn't by rotating around the z-axis. """
    # if first atom doesn't lie in any plane, some extra preparation is required
    atom = group.atoms[atoms_to_put_in_plane[0]]
    
    not_in_any_plane = False
    
    if not atom.x == 0.0 and not atom.y == 0.0 and not atom.z == 0.0:
        not_in_any_plane = True
    
    molecule = rotate_group(group=group, atom=atoms_to_put_in_plane[0], ax="z", not_in_any_plane=not_in_any_plane)
   
    print("Coordinates after first rotation")
    print(molecule)
    # groups_to_plot.append(copy.deepcopy(molecule))
    # labels.append("rotation1")

    """ Second rotation: puts first atom on x-axis by rotating around y-axis. """
    molecule = rotate_group(group=group, atom=atoms_to_put_in_plane[0], ax="y", not_in_any_plane=False)

    print("Coordinates after second rotation")
    print(molecule)
    groups_to_plot.append(copy.deepcopy(molecule))
    labels.append("rotation2")

    """ Third rotation: puts second atom in x-y plane by rotating around the x-axis. """
    molecule = rotate_group(group=group, atom=atoms_to_put_in_plane[1], ax="x", not_in_any_plane=True)

    print("Coordinates after third rotation")
    print(molecule)
    groups_to_plot.append(copy.deepcopy(molecule))
    labels.append("rotation3")

    molecule.invert_if_neccessary()

    print("Coordinates after inversion")
    print(molecule)
    groups_to_plot.append(copy.deepcopy(molecule))
    labels.append("inversion")

    if plot:
        plot_groups(groups_to_plot, labels=labels, bonds=bonds)

    return molecule

def plot_groups(groups, labels, bonds):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    colors = ["red", "gold", "green", "dodgerblue", "fuchsia"]

    for i, group in enumerate(groups):

        atom = list(group.atoms.values())[0]
        ax.scatter(atom.x, atom.y, atom.z, color=colors[i], label=labels[i])
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black') 
        
        for atom in list(group.atoms.values())[1:]:
            ax.scatter(atom.x, atom.y, atom.z, color=colors[i])
            ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black')                 
        
        for bond in bonds:
            x = [group.atoms[bond[0]].x, group.atoms[bond[1]].x]
            y = [group.atoms[bond[0]].y, group.atoms[bond[1]].y]
            z = [group.atoms[bond[0]].z, group.atoms[bond[1]].z]

            ax.plot(x, y, z, color=colors[i])

    ax.legend()
    ax.set_xlabel('X')
    ax.set_xlim(-0.2, 0.2)
    ax.set_ylabel('Y')
    ax.set_ylim(-0.2, 0.2)
    ax.set_zlabel('Z')
    ax.set_zlim(-0.2, 0.2)
    
    plt.show()


if __name__ == "__main__":
    main()


