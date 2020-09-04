from headers import *

from Atom import Atom
from Molecule import Molecule

def load_molecule(filename):
    """ Loads a molecule from a CIF file. Filename should be provided. """  

    with open(filename) as inputfile:
        cif_file = inputfile.readlines()

    reading_coordinates = False
    reading_parameters = False

    target_atoms = []

    molecule = None

    for line in cif_file:

        if line.startswith("_database_code_CSD"):
            label = line.split()[1]
            molecule = Molecule(label)

        # switch reading coordinates
        if line.startswith("_atom_site_label"):
            reading_coordinates = True

        if reading_coordinates == True and line.startswith("loop"):
            reading_coordinates = False

        if reading_coordinates == True and not line.startswith("_"):
            information = line.split()
            atom = Atom(label=information[0], atomtype=information[1], x=information[2], y=information[3], z=information[4])
            molecule.all_atoms.append(atom)

        # switch reading extra parameters
        if line.startswith("_ccdc_geom_distance_query_id"):
            reading_parameters = True

        if reading_parameters and not "END" in line and not line.startswith("_"):
            param = line.split()
            label1 = param[5]
            label2 = param[6]

            if not label1 in target_atoms:
                target_atoms.append(label1)
            if not label2 in target_atoms:
                target_atoms.append(label2)

    # find groups

    # TODO: for now add as one fragment
    molecule.add_fragment(target_atoms)

    return molecule

def count_fragments(molecule, atoms_per_fragment):
    assert (type(len(molecule.highlighted_atoms) / atoms_per_fragment) == int), "amount of fragments should be an integer"
    return len(molecule.highlighted_atoms) / atoms_per_fragment