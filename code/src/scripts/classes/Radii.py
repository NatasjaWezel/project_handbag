import pandas as pd


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

    def get_vdw_distance_contact(self, to_count_contact):
        if to_count_contact == "centroid":
            return self.get_vdw_radius("C")

        # else return vdw radius of the atom the user is interested in
        return self.get_vdw_radius(to_count_contact)

