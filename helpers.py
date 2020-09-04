from headers import *

from collections import defaultdict

from Atom import Atom
from Molecule import Molecule

def load_molecule(filename):
    """ Loads a molecule from a CIF file. Filename should be provided. """  

    with open(filename) as inputfile:
        cif_file = inputfile.readlines()

    reading_coordinates = False
    reading_parameters = False
    reading_bonds = False

    parameterlines = []
    bond_lines = []

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
            parameterlines.append(line)

        if line.startswith("_geom_bond_atom_site_label_1"):
            reading_bonds = True

        if reading_bonds == True and line.startswith("loop"):
            reading_bonds = False

        if reading_bonds == True and not line.startswith("_"):
            bond_lines.append(line)

    molecule = process_parameter_lines(molecule, parameterlines)
    molecule = process_bond_lines(molecule, bond_lines)

    return molecule

def process_bond_lines(molecule, bond_lines):

    for line in bond_lines:
        print(line)
        information = line.split()
        bonded_atom1 = information[0]
        bonded_atom2 = information[1]

        for fragment in molecule.fragments:
            if bonded_atom1 in fragment.atoms.keys():
                fragment.add_bond([bonded_atom1, bonded_atom2])

    return molecule


def process_parameter_lines(molecule, parameterlines):

    target_atoms = defaultdict(list)

    for line in parameterlines:
        param = line.split()
        
        fragment_id = param[1]
        
        label1 = param[5]
        label2 = param[6]

        target_atoms[fragment_id].append(label1)
        target_atoms[fragment_id].append(label2)

    # make sure all labels only appear once
    for fragment_id in target_atoms.keys():
        target_atoms[fragment_id] = list(set(target_atoms[fragment_id]))

    molecule.add_fragments(target_atoms)

    return molecule

def count_fragments(molecule, atoms_per_fragment):
    assert (type(len(molecule.highlighted_atoms) / atoms_per_fragment) == int), "amount of fragments should be an integer"
    return len(molecule.highlighted_atoms) / atoms_per_fragment