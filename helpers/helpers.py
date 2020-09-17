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

        # plot first atom of the fragment and label it
        atom = list(fragment.atoms.values())[0]
        ax.scatter(atom.x, atom.y, atom.z, color=colors[i % len(colors)], label=labels[i])
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black') 
        
        ax = plot_atoms_bonds(i, ax, colors, fragment)

    ax.legend()
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    
    plt.show()

def plot_atoms_bonds(i, ax, colors, fragment):
    """ This function plots the atoms and bonds of a fragment. """ 

    for atom in list(fragment.atoms.values())[1:]:
        ax.scatter(atom.x, atom.y, atom.z, color=colors[i % len(colors)])
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black')                 
    
    for bond in fragment.bonds:
        # TODO: look into this cuz fragments.atoms.keys are supposed to be labels and bond[0] is supposed to be an atom
        if bond[0] in fragment.atoms.keys() and bond[1] in fragment.atoms.keys():
            x = [fragment.atoms[bond[0]].x, fragment.atoms[bond[1]].x]
            y = [fragment.atoms[bond[0]].y, fragment.atoms[bond[1]].y]
            z = [fragment.atoms[bond[0]].z, fragment.atoms[bond[1]].z]

            ax.plot(x, y, z, color=colors[i])
    
    return ax


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

            atom = check_if_label_exists(atom, fragment)

            fragment.add_atom(atom)
            
    fragments.append(fragment)
    
    return fragments

def check_if_label_exists(atom, fragment):
    """ Checks if label already exists in the fragment. Adds the laetter 'a' to it if it does.
        #TODO: should maybe change from a to b, etc. """
        
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
