class Atom:
    """ This class holds atom objects. """ 

    def __init__(self, label, x, y, z):
        self.label = label
        self.distance_to_center = None
        self.part_of = "f"

        if not type(x) == float and "(" in x:
            x = x.split("(")[0]
        if not type(y) == float and "(" in y:
            y = y.split("(")[0]
        if not type(z) == float and "(" in z:
            z = z.split("(")[0]

        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

        self.is_part_of_target = False

    def add_fragment_id(self, fragment_id):
        self.fragment_id = fragment_id

    def __str__(self):
        return self.label + ": " + str(self.x) + ", " + str(self.y) + ", " + str(self.z)

