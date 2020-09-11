from helpers.Fragment import Fragment

import copy
import csv

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
            fragment = Fragment(self.label, fragment_id)

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

  