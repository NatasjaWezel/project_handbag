import time
import numpy as np

from rotation_helpers import perform_rotations

from helpers import load_molecule

def main():

    filenames = ["data/ABOKEJ.NO3_CO_vdw5.cif", "data/AJOWIG.NO3_CO_vdw5.cif", "data/ABINAB.NO3_CO_vdw5.cif"]
    
    results_file_name = "results/saved_data" + str(time.time()) + ".csv"

    for filename in filenames:
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
        molecule.save_fragments_data(results_file_name)





if __name__ == "__main__":
    main()


