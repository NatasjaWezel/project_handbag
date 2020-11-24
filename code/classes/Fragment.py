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

    def center_coordinates(self, center_atom_label):
        """ This is a function that puts any atom you want at the origin of the
            xyz coordinate system, and moves important atoms according to the change. """

        self.center_atom = self.atoms[center_atom_label]

        move_x, move_y, move_z = -self.center_atom.x, -self.center_atom.y, -self.center_atom.z
        self.center_atom.x, self.center_atom.y, self.center_atom.z = 0, 0, 0

        for atom in self.atoms.values():
            if atom is not self.center_atom:
                atom.x += move_x
                atom.y += move_y
                atom.z += move_z

    def __str__(self):
        molecule_string = "Atoms in fragment: " + str(self.id) + "\n"

        for atom in self.atoms.values():
            molecule_string += atom.label + ": " + str(atom.x) + ", " + str(atom.y) + ", " + str(atom.z) + "\n"

        return molecule_string
