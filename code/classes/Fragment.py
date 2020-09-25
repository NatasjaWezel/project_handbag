import math

import pandas as pd

class Fragment:
    def __init__(self, from_entry, fragment_id):
        self.from_entry = from_entry
        self.id = from_entry + str(fragment_id)

        self.atoms = {}


    def add_atom(self, atom):
        self.atoms[atom.label] = atom

    def define_central_group(self, settings):
        self.distances = {}
        atoms = list(self.atoms.values())

        group1 = []
        group2 = []

        i = 0
        # calculate distances
        for atom1 in self.atoms.values():
            
            if len(group1) == 0:
                group1.append(atom1)

            i+=1
            for atom2 in atoms[i:]:
                if not atom1.label + "-" + atom2.label in self.distances.keys():
                    diffx = atom1.x - atom2.x
                    diffy = atom1.y - atom2.y
                    diffz = atom1.z - atom2.z

                    distance = math.sqrt(diffx**2 + diffy**2 + diffz**2)

                    if distance < 1:
                        group1.append(atom2)

                    self.distances[atom1.label + "-" + atom2.label] = distance
        print(self.distances)
                

    def set_center(self, settings):
        for atom in self.atoms.values():
            if atom_type == atom.symbol:
                self.center_atom = atom
                atom.distance_to_center = 0
                atom.part_of = "c"

    def set_vdw_radii(self, filename):
        radii_df = pd.read_csv(filename, header=None)
        radii_df.columns = ["name", "symbol", "radius"]

        for atom in self.atoms.values():
            atom.vdw_radius = float(radii_df[radii_df.symbol == atom.symbol].radius)
    
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
        molecule_string = "Atoms in fragment: " + str(self.id) + "\n"

        for atom in self.atoms.values():
            molecule_string += atom.label + ": " + str(atom.x) + ", " + str(atom.y) + ", " + str(atom.z) + "\n"

        return molecule_string
