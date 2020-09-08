      
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
