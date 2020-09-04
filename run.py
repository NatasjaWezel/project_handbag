from headers import *

import numpy as np

import copy as copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from rotation_helpers import perform_rotations

from tests import check_new_fragment_alignment

from helpers import load_molecule

def main():

    filenames = ["data/NO3_CO_vdw5/ABOKEJ.CO_NO3_vdw5.cif", "data/NO3_CO_vdw5/AJOWIG.CO_NO3_vdw5.cif"]
    bonds = [BONDS_ABOKEJ, BONDS_AJOWIG]

    for bond, filename in zip(bonds, filenames):
        molecule = load_molecule(filename=filename)

        # center on N atom
        atom_to_center = "N1"

        # for fragment in fragments
        fragment = molecule.fragments[0]
        fragment.center_coordinates(atom_to_center=atom_to_center)
        
        # TODO: abstract these from file so you only have to say "center on N atom"
        # now it can give a key error
        atoms_to_put_in_plane = ["O1", "O2"]
        fragment = perform_rotations(fragment, atoms_to_put_in_plane, plot=False, bonds=bond)

        # test if everything went right (if the math is allright)
        print(molecule.label, end=": ")
        check_new_fragment_alignment(fragment, atom_to_center, atoms_to_put_in_plane)


if __name__ == "__main__":
    main()


