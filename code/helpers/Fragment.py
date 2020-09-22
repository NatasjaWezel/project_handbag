import math

import pandas as pd

class Fragment:
    def __init__(self, from_entry, fragment_id):
        self.from_entry = from_entry
        self.fragment_id = fragment_id
        self.unique_id = from_entry + fragment_id

        self.atoms = {}
        self.bonds = []
        self.color = "red"

    def add_atom(self, atom):
        self.atoms[atom.label] = atom

    def add_bond(self, bond):
        self.bonds.append(bond)

    def set_center(self, atom_type):
        for atom in self.atoms.values():
            if atom_type in atom.label:
                self.center_atom = atom
                atom.distance_to_center = 0
                atom.part_of = "c"

    def set_vdw_radii(self, filename):
        radii_df = pd.read_csv(filename, header=None)
        radii_df.columns = ["name", "symbol", "radius"]

        for atom in self.atoms.values():
            atom.vdw_radius = float(radii_df[radii_df.symbol == atom.symbol].radius)

    def find_bonds_NO3_and_distances(self):
        for atom in self.atoms.values():
            
            if not atom == self.center_atom:
                diffx = abs(atom.x - self.center_atom.x)
                diffy = abs(atom.y - self.center_atom.y)
                diffz = abs(atom.z - self.center_atom.z)

                distance = math.sqrt(diffx**2 + diffy**2 + diffz**2)
                atom.distance_to_center = distance

        self.attached_atoms = sorted(self.atoms.values(), key=lambda x: x.distance_to_center, reverse=False)
        
        # TODO: this is now NO3 CO specific. Where to get the bonds from?
        # the first three of these are bonded to the NO3:
        for atom in self.attached_atoms[1:4]:
            self.add_bond([self.center_atom.label, atom.label])
            atom.part_of = "c"

    
    def find_atoms_for_plane(self):
        self.find_bonds_NO3_and_distances()
        return self.attached_atoms[1:3]

    def invert_if_neccessary(self):
        """ Inverts the whole molecule if most of it is on the negative z-axis.
            Does this by mirroring the sign of the z coordinate. """

        z_mean = 0.0

        for atom in self.atoms.values():
            # TODO: make this not-NO3 specific
            # if not "C" in atom.label:
            z_mean += atom.z
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
        molecule_string = "Atoms in fragment: " + str(self.unique_id) + "\n"
        # molecule_string = ""
        for atom in self.atoms.values():
            molecule_string += atom.label + ": " + str(atom.x) + ", " + str(atom.y) + ", " + str(atom.z) + "\n"

        molecule_string += "Bonds in fragment: \n"
        for bond in self.bonds:
            molecule_string += bond[0] + "-" + bond[1] +"\n"


        return molecule_string
