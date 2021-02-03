import pandas as pd
import math


class Radii():
    def __init__(self, RADII_CSV):
        self.radii_filename = RADII_CSV
        self.vdw_radii = {}
        self.cov_radii = {}

    def get_vdw_radius(self, symbol):
        if symbol not in self.vdw_radii.keys():
            radii_df = pd.read_csv(self.radii_filename)

            vdw_radius = float(radii_df[radii_df.symbol == symbol].vdw_radius)

            # TODO: what if it's not in there? like R
            self.vdw_radii[symbol] = vdw_radius

        return self.vdw_radii[symbol]

    def get_cov_radius(self, symbol):
        if symbol not in self.cov_radii.keys():
            radii_df = pd.read_csv(self.radii_filename)

            cov_radius = float(radii_df[radii_df.symbol == symbol].cov_radius)

            # TODO: what if it's not in there? like R
            self.cov_radii[symbol] = cov_radius

        return self.cov_radii[symbol]

    def get_vdw_distance_contact(self, df, settings):
        if settings.to_count_contact == "centroid":
            return self.get_vdw_radius("C")
            # return self.calculate_longest_vdw_radius_contact(df, settings)

        # else return vdw radius of the atom the user is interested in
        return self.get_vdw_radius(settings.to_count_contact)

    def calculate_longest_vdw_radius_contact(self, df, settings):
        # TODO: if there's an R, kick that one out
        longest_distance = 0
        atom_a = None

        # take the first fragment and it's centroid
        first_fragment_df = df[df.fragment_id == df.fragment_id.unique()[0]]
        centroid = first_fragment_df.groupby('fragment_id').mean()

        for _, atom in first_fragment_df.iterrows():
            if atom.label == '-':
                distance = math.sqrt((atom.x - centroid.x)**2 + (atom.y - centroid.y)**2 +
                                     (atom.z - centroid.z)**2)

                if distance > longest_distance:
                    longest_distance = distance
                    atom_a = atom

        longest_vdw_distance = (longest_distance + self.get_vdw_radius(atom_a.symbol))

        return longest_vdw_distance
