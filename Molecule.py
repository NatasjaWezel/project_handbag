from headers import *

# needed to plot the molecule (in 3D)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# needed for matrix calculations
import numpy as np

import copy

class Molecule:
    """ This class can hold a molecule, which is a list of atoms. """  

    def __init__(self, label):
        self.label = label
        
        self.all_atoms = []
        self.fragments = []

    def extend_molecule(self, atom):
        self.all_atoms.append(atom)
    
    def add_fragments(self, fragment_labels):
        """ This function adds a pair of target and coordination group. """ 
        
        for fragment_id in fragment_labels.keys():
            fragment = Fragment(fragment_id)

            for atom in self.all_atoms:
                if atom.label in fragment_labels[fragment_id]:
                    # this can't be a reference since atoms can be in multiple
                    # fragments and we're gonna manipulate its coordinates
                    fragment.atoms[atom.label] = copy.deepcopy(atom)

            self.fragments.append(fragment)

    def center_fragments(self, atom_to_center):
        for fragment in self.fragments:
            fragment.set_center(atom_to_center)
            fragment.center_coordinates()

    def __str__(self):
        molecule_string = "\nFragments in molecule: " + self.label + "\n"

        for fragment in self.fragments:
            molecule_string += "\n" + str(fragment)

        return molecule_string

        
class Fragment:
    def __init__(self, fragment_id):
        self.fragment_id = fragment_id
        self.atoms = {}
        self.bonds = []

    def add_bond(self, bond):
        self.bonds.append(bond)

    def set_center(self, atom_type):
        for atom in self.atoms.values():
            if atom_type in atom.label:
                self.center_atom = atom

    def find_atoms_for_plane(self):
        attached_atoms = []

        for bond in self.bonds:
            if bond[0] == self.center_atom.label:
                attached_atoms.append(bond[1])
            elif bond[1] == self.center_atom.label:
                attached_atoms.append(bond[0])

        return attached_atoms[0:2]

    def invert_if_neccessary(self):
        """ Inverts the whole molecule if most of it is on the negative z-axis.
            Does this by mirroring the sign of the z coordinate. """

        z_mean = 0.0

        for atom in self.atoms.values():
            z_mean += atom.x
        z_mean = z_mean/len(self.atoms)

        # switch signs
        if z_mean < 0:
            for atom in self.atoms.values():
                if atom.z < 0:
                    atom.z = abs(atom.z)
                elif atom.z > 0:
                    atom.z = atom.z - 2*atom.z

    def center_coordinates(self):
        """ This is a function that puts any atom you want at the origin of the
            xyz coordinate system, and moves important atoms according to the change. """ 

        move_x, move_y, move_z = -self.center_atom.x, -self.center_atom.y, -self.center_atom.z 
        self.center_atom.x, self.center_atom.y, self.center_atom.z = 0, 0, 0

        for atom in self.atoms.values():
            if not atom is self.center_atom:
                atom.x += move_x
                atom.y += move_y
                atom.z += move_z

    def __str__(self):
        molecule_string = "Atoms in fragment: " + self.fragment_id + "\n"

        for atom in self.atoms.values():
            x = 0.0 if atom.x < CUT_OFF_ZERO and atom.x > -CUT_OFF_ZERO else atom.x 
            y = 0.0 if atom.y < CUT_OFF_ZERO and atom.y > -CUT_OFF_ZERO else atom.y 
            z = 0.0 if atom.z < CUT_OFF_ZERO and atom.z > -CUT_OFF_ZERO else atom.z 

            molecule_string += atom.label + ": " + str(x) + ", " + str(y) + ", " + str(z) + "\n"

        molecule_string += "Bonds in fragment: \n"
        for bond in self.bonds:
            molecule_string += bond[0] + "-" + bond[1] +"\n"

        return molecule_string
