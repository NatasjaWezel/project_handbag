from Atom import Atom
from Molecule import Molecule

from rotation_helpers import rotate_molecule, find_angles

def main():
    molecule = Molecule()

    atom1 = Atom("C", "C", 0.15, 0.06, 0.04)
    atom2 = Atom("O", "O", 0.0, 0.0, 0.0)

    molecule.extend_molecule(atom1)
    molecule.extend_molecule(atom2)

    molecule.highlight_target_atoms(labels=["O", "C"])

    molecule.plot_molecule(BONDS=[["O", "C"]])

    print("\nOriginal coordinates")
    for atom in molecule.highlighted_atoms.values():
        print(atom)
    print()

    # project atom onto xy plane and find angle with x-axis
    alpha, _, _ = find_angles(coord_vector=[atom1.x, atom1.y, 0.0])
    # rotate this angle around the z-axis in the positive direction, to put the atom in the xz plane.
    molecule = rotate_molecule(molecule, angle=alpha, ax="z")

    print("After first rotation")
    for atom in molecule.highlighted_atoms.values():
        print(atom)
    print()

    molecule.plot_molecule(BONDS=[["O", "C"]])

    # find the angle with the y-axis
    alpha, _, _ = find_angles(coord_vector=[atom1.x, atom1.y, atom1.z])
    # rotate this angle around the y-axis to put the atom on the y axis.
    molecule = rotate_molecule(molecule, angle=alpha, ax="y")

    print("After second rotation")
    for atom in molecule.highlighted_atoms.values():
        print(atom)
    print()
    
    molecule.plot_molecule(BONDS=[["O", "C"]])




if __name__ == "__main__":
    main()