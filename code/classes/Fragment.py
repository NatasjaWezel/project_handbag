import math

import pandas as pd

class Fragment:
    def __init__(self, from_entry, fragment_id):
        self.from_entry = from_entry
        self.id = from_entry + str(fragment_id)
        
        self.bonds = []
        self.color = "red"
        self.atoms = {}

    def add_atom(self, atom):
        self.atoms[atom.label] = atom

    def define_central_group(self, atom):
        atom.add_to_central_group()

    def find_atoms_for_plane(self, settings):
        plane_atoms = []

        for key, value in settings.central_group_atoms.items():
            if value == 1 and settings.center_ring:
                for atom in self.atoms.values():
                    if atom.in_central_group and atom is not self.center_atom and atom.symbol == key:
                        plane_atoms.append(atom)

        for atom in self.atoms.values():
            if atom.in_central_group and atom is not self.center_atom and not atom in plane_atoms:
                plane_atoms.append(atom)
        
        return plane_atoms[:2]

    def invert_if_neccessary(self):
        """ Inverts the whole molecule if most of it is on the negative z-axis.
            Does this by mirroring the sign of the z coordinate. """

        z_mean = 0.0

        for atom in self.atoms.values():
            z_mean += atom.z
        z_mean = z_mean/len(self.atoms)

        # switch signs
        if z_mean < 0:
            for atom in self.atoms.values():
                if atom.z < 0:
                    atom.z = abs(atom.z)
                elif atom.z > 0:
                    atom.z = atom.z - 2*atom.z

    def find_moves_xyz_to_center(self, center_ring):
        if center_ring != []:
            no_atoms = len(center_ring)
            move_x = -sum([atom.x for atom in center_ring])/no_atoms
            move_y = -sum([atom.y for atom in center_ring])/no_atoms
            move_z = -sum([atom.z for atom in center_ring])/no_atoms
        elif self.center_atom:
            move_x, move_y, move_z = -self.center_atom.x, -self.center_atom.y, -self.center_atom.z 

            self.center_atom.x, self.center_atom.y, self.center_atom.z = 0, 0, 0
        else:
            assert center_ring != [] or self.center_atom, "No center atom and no center of ring found.... :/"

        return move_x, move_y, move_z

    def center_coordinates(self, settings):
        """ This is a function that puts any atom you want at the origin of the
            xyz coordinate system, and moves important atoms according to the change. """ 
        
        self.center_atom = None
        center_ring = []

        for atom in self.atoms.values():

            if self.center_atom == None and atom.in_central_group and atom.symbol == settings.center_atom:
                self.center_atom = atom
            elif atom.in_central_group and atom.symbol == settings.center_ring:
                center_ring.append(atom)
        
        move_x, move_y, move_z = self.find_moves_xyz_to_center(center_ring)

        for atom in self.atoms.values():
            if not atom is self.center_atom:
                atom.x += move_x
                atom.y += move_y
                atom.z += move_z

    def __str__(self):
        molecule_string = "Atoms in fragment: " + str(self.id) + "\n"

        for atom in self.atoms.values():
            molecule_string += atom.label + ": " + str(atom.x) + ", " + str(atom.y) + ", " + str(atom.z) + "\n"

        return molecule_string
