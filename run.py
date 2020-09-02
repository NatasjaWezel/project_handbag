from Atom import Atom
from Molecule import Molecule

import numpy as np

import copy as copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from rotation_helpers import rotate_molecule

BONDS_ABOKEJ = [["N1", "O1"], 
                ["N1", "O2"], 
                ["N1", "O2A"],
                ["C5", "O6"]]

def main():

    molecule = load_molecule(filename="data/NO3_CO_vdw5/ABOKEJ.CO_NO3_vdw5.cif")

    # abstract these from file so you only have to say "center on N atom"
    atom_to_center = "N1"
    atoms_to_put_in_plane = ["O2", "O2A"]

    molecule.center_coordinates(atom_to_center=atom_to_center)
    
    molecule = perform_rotations(molecule, atoms_to_put_in_plane, plot=False)


def perform_rotations(molecule, atoms_to_put_in_plane, plot):
    print("Original coordinates")
    print(molecule)
    
    molecules_to_plot = []
    molecules_to_plot.append(copy.deepcopy(molecule))

    molecule = rotate_molecule(molecule=molecule, atom=atoms_to_put_in_plane[0], ax="z", not_in_any_plane=True)

    print("Coordinates after first rotation")
    print(molecule)
    molecules_to_plot.append(copy.deepcopy(molecule))

    molecule = rotate_molecule(molecule=molecule, atom=atoms_to_put_in_plane[0], ax="y", not_in_any_plane=False)

    print("Coordinates after second rotation")
    print(molecule)
    molecules_to_plot.append(copy.deepcopy(molecule))

    molecule = rotate_molecule(molecule=molecule, atom=atoms_to_put_in_plane[1], ax="x", not_in_any_plane=False)

    print("Coordinates after third rotation")
    print(molecule)
    molecules_to_plot.append(copy.deepcopy(molecule))

    if plot:
        plot_molecules([molecule])

    return

def plot_molecules(molecules):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    colors = ["red", "green", "blue", "purple"]

    for i, molecule in enumerate(molecules):
        for atom in molecule.highlighted_atoms.values():
            ax.scatter(atom.x, atom.y, atom.z, color=colors[i])
            ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black')                 
        for bond in BONDS_ABOKEJ:
            x = [molecule.highlighted_atoms[bond[0]].x, molecule.highlighted_atoms[bond[1]].x]
            y = [molecule.highlighted_atoms[bond[0]].y, molecule.highlighted_atoms[bond[1]].y]
            z = [molecule.highlighted_atoms[bond[0]].z, molecule.highlighted_atoms[bond[1]].z]

            ax.plot(x, y, z, color=colors[i])


    ax.set_xlabel('X')
    ax.set_xlim(-0.2, 0.2)
    ax.set_ylabel('Y')
    ax.set_ylim(-0.2, 0.2)
    ax.set_zlabel('Z')
    ax.set_zlim(-0.2, 0.2)
    
    plt.show()


def load_molecule(filename):
    """ Loads a molecule from a CIF file. Filename should be provided. """  

    with open(filename) as inputfile:
        cif_file = inputfile.readlines()

    reading_coordinates = False
    reading_parameters = False
    target_atoms = []

    molecule = Molecule()

    for line in cif_file:

        # switch reading coordinates
        if line.startswith("_atom_site_label"):
            reading_coordinates = True

        if reading_coordinates == True and line.startswith("loop"):
            reading_coordinates = False

        if reading_coordinates == True and not line.startswith("_"):
            information = line.split()
            atom = Atom(label=information[0], atomtype=information[1], x=information[2], y=information[3], z=information[4])
            molecule.extend_molecule(atom)

        # switch reading extra parameters
        if line.startswith("_ccdc_geom_distance_query_id"):
            reading_parameters = True

        if reading_parameters and not "END" in line and not line.startswith("_"):
            param = line.split()
            label1 = param[5]
            label2 = param[6]

            if not label1 in target_atoms:
                target_atoms.append(label1)
            if not label2 in target_atoms:
                target_atoms.append(label2)

    molecule.highlight_target_atoms(target_atoms)
    
    return molecule





if __name__ == "__main__":
    main()


