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

    def define_central_group(self, settings):
        distances, group1 = self.calculate_distances(settings)

        central_group = find_central_group(distances, list(self.atoms.values()), group1, settings)

        if not central_group == None:
            for atom in central_group:
                atom.add_to_central_group()
        
        return central_group    

    def calculate_distances(self, settings):
        distances = {}
        atoms = list(self.atoms.values())

        group1 = []

        i = 0

        for atom1 in atoms:
            if group1 == []:
                group1.append(atom1)
            
            i += 1
            for atom2 in atoms[i:]:
                diffx = atom1.x - atom2.x
                diffy = atom1.y - atom2.y
                diffz = atom1.z - atom2.z

                distance = math.sqrt(diffx**2 + diffy**2 + diffz**2)

                distances[atom1.label + "-" + atom2.label] = distance

                # TODO: this won't always work for phenyl for exmaple. Make something recursive
                # print(atom1.label + "-" + atom2.label, distance, settings.get_cov_radius(atom1.symbol) + settings.get_cov_radius(atom2.symbol))
                if distance < settings.get_cov_radius(atom1.symbol) + settings.get_cov_radius(atom2.symbol) - 0.01:
                    # print("BOND!")
                    if atom1 in group1 and not atom2 in group1:
                        group1.append(atom2)
                    elif atom2 in group1 and not atom1 in group1:
                        group1.append(atom1)
        
        return distances, group1

    def find_atoms_for_plane(self):
        plane_atoms = []
        for atom in self.atoms.values():
            if atom.in_central_group and atom is not self.center_atom:
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
            if atom.in_central_group and atom.symbol == settings.center_atom:
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

def count_atoms(group):
    atoms_count = {}
    for atom in group:
        if atom.symbol in atoms_count.keys():
            atoms_count[atom.symbol] += 1
        else:
            atoms_count[atom.symbol] = 1
    return atoms_count

# TODO: find out why there are extra atoms added sometimes.
def find_central_group(distances, atoms, group1, settings, tolerance=0.0):
    
    group2 = []
    for atom in atoms:
        if atom not in group1:
            group2.append(atom)

    atom_count1 = count_atoms(group1)
    atom_count2 = count_atoms(group2)

    if settings.central_group_atoms == atom_count1:
        return group1
    elif settings.central_group_atoms == atom_count2:
        return group2
    else:
        if tolerance >= 1:
            # print("Tolerance over 1. Something must've gone wrong")
           
            return None
 
        group1 = new_group1(distances, atoms, settings, tolerance)
        tolerance += 0.01
        return find_central_group(distances, atoms, group1, settings, tolerance)

def new_group1(distances, atoms, settings, tolerance):
    group1 = []
    i = 0

    for atom1 in atoms:
        if group1 == []:
            group1.append(atom1)
        
        i += 1
        for atom2 in atoms[i:]:

                distance = distances[atom1.label + "-" + atom2.label]
                
                # TODO: this won't always work for phenyl for exmaple. Make something recursive
                if distance < settings.get_cov_radius(atom1.symbol) + settings.get_cov_radius(atom2.symbol) + tolerance:
                    if atom1 in group1 and not atom2 in group1:
                        group1.append(atom2)
                    elif atom2 in group1 and not atom1 in group1:
                        group1.append(atom1)
    return group1
