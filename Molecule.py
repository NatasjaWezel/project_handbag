from headers import *

# needed to plot the molecule (in 3D)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# needed for matrix calculations
import numpy as np

class Molecule:
    """ This class can hold a molecule, which is a list of atoms. """  

    def __init__(self, label):
        self.label = label
        
        self.all_atoms = []
        self.groups = []

    def extend_molecule(self, atom):
        self.all_atoms.append(atom)
    
    def add_group(self, group_labels):
        """ This function adds a pair of target and coordination group. """ 
        group = Group()

        for atom in self.all_atoms:
            if atom.label in group_labels:
                group.atoms[atom.label] = atom

        self.groups.append(group)


class Group:
    def __init__(self):
        self.atoms = {}

    def invert_if_neccessary(self):
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

    def center_coordinates(self, atom_to_center):
        """ This is a function that puts any atom you want at the origin of the
            xyz coordinate system, and moves important atoms according to the change. """ 

        atom_to_center = self.atoms[atom_to_center]

        move_x, move_y, move_z = -atom_to_center.x, -atom_to_center.y, -atom_to_center.z 
        atom_to_center.x, atom_to_center.y, atom_to_center.z = 0, 0, 0

        for atom in self.atoms.values():
            if not atom is atom_to_center:
                atom.x += move_x
                atom.y += move_y
                atom.z += move_z

    def __str__(self):
        molecule_string = ""

        for atom in self.atoms.values():
            x = 0.0 if atom.x < CUT_OFF_ZERO and atom.x > -CUT_OFF_ZERO else atom.x 
            y = 0.0 if atom.y < CUT_OFF_ZERO and atom.y > -CUT_OFF_ZERO else atom.y 
            z = 0.0 if atom.z < CUT_OFF_ZERO and atom.z > -CUT_OFF_ZERO else atom.z 

            molecule_string += atom.label + ": " + str(x) + ", " + str(y) + ", " + str(z) + "\n"

        return molecule_string
