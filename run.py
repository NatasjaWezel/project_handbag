import time
import numpy as np

from rotation_helpers import perform_rotations

from helpers import load_molecule

import pickle

import sys

def main():
    
    if len(sys.argv) != 3:
        print("Usage: python run.py <inputfilename> <outputfilename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    outputfilename = sys.argv[2]
    
    #TODO: think of a smarter way than loading/saving the class for every entry?
    # like a csv file, but that missed the bonds?
    try:
        f = open(outputfilename, 'rb')
        saved_molecules = pickle.load(f)
        f.close()
    except FileNotFoundError:
        saved_molecules = Save_molecules()

    molecule = load_molecule(filename=filename)

    # center on N atom
    atom_to_center = "N"
    molecule.center_fragments(atom_to_center)

    for fragment in molecule.fragments:
        atoms_to_put_in_plane = fragment.find_atoms_for_plane()
        
        fragment = perform_rotations(fragment, atoms_to_put_in_plane, plot=False)
        print(molecule.label, fragment.fragment_id, "Passed all checks. Rotation OK")
        
        fragment.invert_if_neccessary()

    # per file save the new datapoints of each fragment
    saved_molecules.add_molecule(molecule)

    save_file = open(outputfilename, "wb")
    pickle.dump(saved_molecules, save_file)
    save_file.close()


class Save_molecules():

    def __init__(self):
        self.molecules = []

    def add_molecule(self, molecule):
        self.molecules.append(molecule)



if __name__ == "__main__":
    main()


