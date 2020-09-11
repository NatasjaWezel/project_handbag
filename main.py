import time
import numpy as np

from helpers.rotation_helpers import perform_rotations

from helpers.helpers import load_molecule

import pickle

import sys

def main():
    
    if len(sys.argv) != 3:
        print("Usage: python run.py <inputfilename> <outputfilename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    outputfilename = sys.argv[2]

    molecule = load_molecule(filename=filename)

    # center on N atom
    atom_to_center = "N"
    molecule.center_fragments(atom_to_center)

    for fragment in molecule.fragments:
        atoms_to_put_in_plane = fragment.find_atoms_for_plane()
        
        print(molecule.label, end=" ")
        fragment = perform_rotations(fragment, atoms_to_put_in_plane, plot=False)
        print(fragment.fragment_id, "Passed all checks. Rotation OK")
        
        fragment.invert_if_neccessary()
        
    molecule.save_fragments_data(filename=outputfilename)
    

if __name__ == "__main__":
    main()


