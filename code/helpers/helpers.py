from collections import defaultdict

from helpers.Atom import Atom
from helpers.Fragment import Fragment

def load_fragments_from_coords(filename):
    """ Loads the list of aligned fragments from a .cor file. """

    with open(filename) as inputfile:
        lines = inputfile.readlines()

    fragments = []
    fragment = None

    for line in lines:
        if "FRAG" in line:
            if fragment:
                fragments.append(fragment)

            information = line.split('**')
            fragment = Fragment(fragment_id=information[2].strip(), from_entry=information[0].strip())
        else:
            information = line.split()
            x, y, z = information[0].split("("), information[1].split("("), information[2].split("(")

            atom = Atom(label=information[0].strip("%"), coordinates=[float(x[0]), float(y[0]), float(z[0])])

            atom = check_if_label_exists(atom, fragment)

            fragment.add_atom(atom)
            
    fragments.append(fragment)
    
    return fragments


def check_if_label_exists(atom, fragment):
    """ Checks if label already exists in the fragment. Adds the laetter 'a' to it if it does. """
        
    if atom.label in fragment.atoms.keys():
        atom.label += 'a'
        atom = check_if_label_exists(atom, fragment)
    
    return atom


def process_bond_lines(molecule, bond_lines):
    """ Adds bonds to each fragment in a molecule, by reading which bonds 
        exist in the inputfile. """

    for line in bond_lines:
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
