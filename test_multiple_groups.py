from helpers import load_molecule

def main():
    filenames = ["data/NO3_CO_vdw5/ABINAB.CO_NO3_vdw5.cif"]

    for bond, filename in zip(bonds, filenames):
        molecule = load_molecule(filename=filename)

        # center on N atom
        atom_to_center = "N1"

        # for group in groups
        # group = molecule.groups[0]
        # group.center_coordinates(atom_to_center=atom_to_center)
        
        # TODO: abstract these from file so you only have to say "center on N atom"
        # now it can give a key error
        # atoms_to_put_in_plane = ["O1", "O2"]
        # group = perform_rotations(group, atoms_to_put_in_plane, plot=False, bonds=bond)

        # test if everything went right (if the math is allright)
        # print(molecule.label, end=": ")
        # check_new_group_alignment(group, atom_to_center, atoms_to_put_in_plane)


if __name__ == "__main__":
    main()