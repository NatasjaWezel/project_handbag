import pandas as pd

class Settings():
    def __init__(self, inputfilename):

        self.outputfilename = "./results/" + inputfilename.rsplit('\\')[-1].rsplit('.', 1)[0] + "_aligned.csv"
        self.working_directory = "X"

    def set_central_group(self, groupname):
        self.central_group_name = groupname

        # TODO make this variable
        groupnames = pd.read_csv("./data/central_groups.txt", header=0)

        df = groupnames[groupnames.name == groupname]

        self.central_group_atoms = {}

        for _, row in df.iterrows():
            self.central_group_atoms[row["atom"]] = int(row["amount"])

            if row.center == "True" and row.amount == "1":
                self.center_atom = row.atom
                self.center_ring = False
            elif row.center == "True":
                self.center_atom = False
                self.center_ring = row.atom
                
        print(self.central_group_atoms)

        