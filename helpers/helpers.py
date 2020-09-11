from collections import defaultdict

from helpers.Atom import Atom
from helpers.Fragment import Fragment
from helpers.Molecule import Molecule

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_fragments(fragments, labels):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    colors = ["red", "sandybrown", "gold", "chartreuse", "green", 
                "mediumturquoise", "dodgerblue", "darkblue", "slateblue",
                "mediumorchid", "fuchsia"]

    for i, fragment in enumerate(fragments):

        atom = list(fragment.atoms.values())[0]
        ax.scatter(atom.x, atom.y, atom.z, color=colors[i % len(colors)], label=labels[i])
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black') 
        
        for atom in list(fragment.atoms.values())[1:]:
            ax.scatter(atom.x, atom.y, atom.z, color=colors[i % len(colors)])
            ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black')                 
        
        for bond in fragment.bonds:
            if bond[0] in fragment.atoms.keys() and bond[1] in fragment.atoms.keys():
                x = [fragment.atoms[bond[0]].x, fragment.atoms[bond[1]].x]
                y = [fragment.atoms[bond[0]].y, fragment.atoms[bond[1]].y]
                z = [fragment.atoms[bond[0]].z, fragment.atoms[bond[1]].z]

                ax.plot(x, y, z, color=colors[i])

    ax.legend()
    ax.set_xlabel('X')
    ax.set_xlim(-2, 6)
    ax.set_ylabel('Y')
    ax.set_ylim(-2, 6)
    ax.set_zlabel('Z')
    ax.set_zlim(-2, 6)
    
    plt.show()


def load_fragments_from_coords(filename):
    with open(filename) as inputfile:
        lines = inputfile.readlines()

    fragments = []
    fragment = None

    for line in lines:
        if "FRAG" in line:
            if fragment:
                fragments.append(fragment)

            information = line.split('**')
            entry = information[0].strip()
            fragment_id = information[2].strip()
            fragment = Fragment(fragment_id=fragment_id, from_entry=entry)
        else:
            information = line.split()
            atom = Atom(label=information[0].strip("%"), x=information[1], y=information[2], z=information[3])
            
            # TODO: something nice recursive here
            if atom.label not in fragment.atoms.keys():
                fragment.add_atom(atom)
            else:
                atom.label += 'a'
                if atom.label not in fragment.atoms.keys():
                    fragment.add_atom(atom)
                else:
                    atom.label += 'b'
                    fragment.add_atom(atom)

    fragments.append(fragment)
    
    return fragments





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
            atom = Atom(label=information[0], x=information[2], y=information[3], z=information[4])
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
