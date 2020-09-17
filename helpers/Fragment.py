from helpers.headers import CUT_OFF_ZERO
import math

class Fragment:
    def __init__(self, from_entry, fragment_id):
        self.from_entry = from_entry
        self.fragment_id = fragment_id
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

        # the other two are also attached
        self.add_bond([self.attached_atoms[4].label, self.attached_atoms[5].label])

    
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
        molecule_string = "Atoms in fragment: " + str(self.fragment_id) + "\n"
        # molecule_string = ""
        for atom in self.atoms.values():
            x = 0.0 if atom.x < CUT_OFF_ZERO and atom.x > -CUT_OFF_ZERO else atom.x 
            y = 0.0 if atom.y < CUT_OFF_ZERO and atom.y > -CUT_OFF_ZERO else atom.y 
            z = 0.0 if atom.z < CUT_OFF_ZERO and atom.z > -CUT_OFF_ZERO else atom.z 

            molecule_string += atom.label + ": " + str(x) + ", " + str(y) + ", " + str(z) + "\n"

        molecule_string += "Bonds in fragment: \n"
        if self.bonds:
            for bond in self.bonds:
                molecule_string += bond[0] + "-" + bond[1] +"\n"
        else:
            molecule_string += "None\n"

        return molecule_string
