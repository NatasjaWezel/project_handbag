from helpers import load_molecule
from rotation_helpers import perform_rotations
from tests import check_new_fragment_alignment

def main():
    filenames = ["data/ABINAB.NO3_CO_vdw5.cif"]

    for filename in filenames:
        molecule = load_molecule(filename=filename)

        # center on N atom
        atom_to_center = "N"

        molecule.center_fragments(atom_to_center)

        for fragment in molecule.fragments:
            atoms_to_put_in_plane = fragment.find_atoms_for_plane()
            
            fragment = perform_rotations(fragment, atoms_to_put_in_plane, plot=False)
            
            print(molecule.label, "fragment", fragment.fragment_id, end=": ")
            check_new_fragment_alignment(fragment, fragment.center_atom.label, atoms_to_put_in_plane)

            fragment.invert_if_neccessary()

if __name__ == "__main__":
    main()