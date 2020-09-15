from helpers.rotation_helpers import perform_rotations
from helpers.helpers import load_fragments_from_coords

import csv
import sys

def main():
    
    if len(sys.argv) != 3:
        print("Usage: python load_from_coords.py <inputfilename> <outputfilename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    outputfilename = sys.argv[2]

    fragments = load_fragments_from_coords(filename=filename)

    rotated_fragments = []

    # center on N atom
    atom_to_center = "N"

    for fragment in fragments:
        try:
            fragment.set_center(atom_to_center)
            fragment.center_coordinates()
            atoms_to_put_in_plane = fragment.find_atoms_for_plane()

            # print(fragment.from_entry, end=" ")
            fragment = perform_rotations(fragment, atoms_to_put_in_plane, plot=False)
            # print(fragment.fragment_id, "Passed all checks. Rotation OK")
            
            fragment.invert_if_neccessary()

            rotated_fragments.append(fragment)
        except AssertionError as msg:
            print(msg)
        
    print(len(rotated_fragments), "/", len(fragments), ' rotated succesfully')

    with open(outputfilename, "a", newline="") as outputfile:
        writer = csv.writer(outputfile)

        for fragment in rotated_fragments:
            for atom in fragment.atoms.values():
                writer.writerow([fragment.from_entry, fragment.fragment_id, atom.label, atom.part_of, atom.x, atom.y, atom.z])
    

if __name__ == "__main__":
    main()


