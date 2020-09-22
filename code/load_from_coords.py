# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the fragments exported from a conquest query and 
# aligns the central groups by using rotation matrices and other linear algebra.
# It then saves the new coordinates in a .csv file.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.rotation_helpers import perform_rotations
from helpers.helpers import load_fragments_from_coords

import csv
import sys

def main():
    
    if len(sys.argv) != 2:
        print("Usage: python load_from_coords.py <path/to/inputfile>")
        sys.exit(1)
    
    filename = sys.argv[1]
    outputfilename = "results/" + filename.rsplit('\\')[-1].rsplit('.', 1)[0] + "_aligned.csv"

    fragments = load_fragments_from_coords(filename=filename)
    
    rotated_fragments = []

    # center on N atom
    atom_to_center = "N"

    for fragment in fragments:
        try:
            fragment.set_center(atom_to_center)        
            fragment.center_coordinates()

            atoms_to_put_in_plane = fragment.find_atoms_for_plane()

            fragment = perform_rotations(fragment, atoms_to_put_in_plane)
            
            fragment.invert_if_neccessary()

            rotated_fragments.append(fragment)
        except AssertionError as msg:
            print(msg)
                    
    print(len(rotated_fragments), "/", len(fragments), ' rotated succesfully')

    with open(outputfilename, "a", newline="") as outputfile:
        writer = csv.writer(outputfile)

        for fragment in rotated_fragments:
            for atom in fragment.atoms.values():
                writer.writerow([fragment.from_entry, fragment.fragment_id, fragment.from_entry + fragment.fragment_id, atom.label, atom.symbol, atom.part_of, atom.x, atom.y, atom.z])
    

if __name__ == "__main__":
    main()



