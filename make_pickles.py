import time
import numpy as np

from helpers.rotation_helpers import perform_rotations
from helpers.helpers import load_molecule

import pickle

import sys

def main():
    outputfilename = "results/saved_data_NO3_CO" + str(time.time()) + ".pkl"

    try:
        f = open(outputfilename, 'rb')
        saved_molecules = pickle.load(f)
        f.close()
    except FileNotFoundError:
        saved_molecules = Save_molecules()

    
    filenames = ["data/NO3_CO_vdw5/ABOKEJ.NO3_CO_vdw5.cif", "data/NO3_CO_vdw5/AJOWIG.NO3_CO_vdw5.cif", "data/NO3_CO_vdw5/ABINAB.NO3_CO_vdw5.cif"]

    for filename in filenames:
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


