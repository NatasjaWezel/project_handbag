import pandas as pd


class Fingerprint():
    def __init__(self, fingerprintcsv, central):
        self.central_group_name = central
        self.csv = fingerprintcsv

        df = pd.read_csv(fingerprintcsv)
        print(df)

    def get_description(self, index):
        pass
