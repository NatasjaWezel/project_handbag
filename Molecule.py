
# needed to plot the molecule (in 3D)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# needed for matrix calculations
import numpy as np


BONDS_ABOKEJ = [["N1", "O1"], 
                ["N1", "O2"], 
                ["N1", "O2A"],
                ["C5", "O6"]]


class Molecule:
    """ This class can hold a molecule, which is a list of atoms. """  

    def __init__(self):
        self.atoms = []
        self.highlighted_atoms = {}

    def extend_molecule(self, atom):
        self.atoms.append(atom)
    
    def highlight_target_atoms(self, labels):
        """ This function highlights the target and coordination group. """ 

        for atom in self.atoms:
            if atom.label in labels:
                atom.highlight()
                self.highlighted_atoms[atom.label] = atom

    
    def center_coordinates(self, atom_to_center):
        """ This is a function that puts any atom you want at the origin of the
            xyz coordinate system, and moves important atoms according to the change. """ 

        atom_to_center = self.highlighted_atoms[atom_to_center]

        move_x, move_y, move_z = -atom_to_center.x, -atom_to_center.y, -atom_to_center.z 
        atom_to_center.x, atom_to_center.y, atom_to_center.z = 0, 0, 0

        for atom in self.highlighted_atoms.values():
            if not atom is atom_to_center:
                atom.x += move_x
                atom.y += move_y
                atom.z += move_z

    def plot_molecule(self, only_highlighted=True, BONDS=False):
        """ Plots the molecule. Uses colors for Nitrogen, oxygen and carbon. In principle
            only plots the important atoms. """ 
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for atom in self.atoms:
            if atom.is_part_of_target:
                color = "grey"
                if "O" in atom.label:
                    color = "red"
                elif "N" in atom.label:
                    color = "blue"
                elif "C" in atom.label:
                    color = "cyan"
                
                ax.scatter(atom.x, atom.y, atom.z, color=color)

            elif only_highlighted == False:
                ax.scatter(atom.x, atom.y, atom.z, color="grey")
        
        if BONDS:
            for bond in BONDS:
                x = [self.highlighted_atoms[bond[0]].x, self.highlighted_atoms[bond[1]].x]
                y = [self.highlighted_atoms[bond[0]].y, self.highlighted_atoms[bond[1]].y]
                z = [self.highlighted_atoms[bond[0]].z, self.highlighted_atoms[bond[1]].z]

            ax.plot(x, y, z, color="grey")


        ax.set_xlabel('X')
        ax.set_xlim(-0.2, 0.2)
        ax.set_ylabel('Y')
        ax.set_ylim(-0.2, 0.2)
        ax.set_zlabel('Z')
        ax.set_zlim(-0.02, 0.2)
        
        plt.show()

    def __str__(self):
        molecule_string = ""

        for atom in self.highlighted_atoms.values():
            molecule_string += atom.label + ": " + str(atom.x) + ", " + str(atom.y) + ", " + str(atom.z) + "\n"

        return molecule_string