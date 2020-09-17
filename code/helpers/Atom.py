class Atom:
    """ This class holds atom objects. """ 

    def __init__(self, label, coordinates):
        self.label = label
        self.distance_to_center = None
        self.part_of = "f"

        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = coordinates[2]

        self.is_part_of_target = False

    def add_fragment_id(self, fragment_id):
        self.fragment_id = fragment_id

    def __str__(self):
        return self.label + ": " + str(self.x) + ", " + str(self.y) + ", " + str(self.z)

