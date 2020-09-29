import pandas as pd

from helpers.headers import RADII_CSV, DATADIR, RESULTSDIR
import os

class Settings():
    def __init__(self, inputfilename):
        
        title = inputfilename.rsplit('\\')[-1].rsplit('.', 1)[0].rsplit('_aligned', 1)[0]
            
        self.output_folder = RESULTSDIR + title + "\\"

        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

        self.outputfile_prefix = self.output_folder + title

        self.radii_filename = RADII_CSV

        self.vdw_radii = {}
        self.cov_radii = {}

    def get_aligned_csv_filename(self):
        aligned_csv_name = self.outputfile_prefix + "_aligned.csv"
        return aligned_csv_name
        
    def get_avg_fragment_filename(self):
        avg_fragment_filename = self.outputfile_prefix + "_avg_fragment.pkl"
        return avg_fragment_filename

    def get_density_df_filename(self, resolution):
        density_df_filename = self.outputfile_prefix + "_" + str(resolution) + "_density.hdf"
        return density_df_filename

    def get_density_plotname(self, resolution):
        density_plotname = self.outputfile_prefix + "_" + str(resolution) + "_density.pdf"
        return density_plotname

    def get_vdw_radius(self, symbol):
        
        if symbol in self.vdw_radii.keys():
            return self.vdw_radii[symbol]
        else:
            radii_df = pd.read_csv(self.radii_filename)

            vdw_radius = float(radii_df[radii_df.symbol == symbol].vdw_radius)

            self.vdw_radii[symbol] = vdw_radius
            return vdw_radius

    def get_cov_radius(self, symbol):
        
        if symbol in self.cov_radii.keys():
            return self.cov_radii[symbol]
        else:
            radii_df = pd.read_csv(self.radii_filename)

            cov_radius = float(radii_df[radii_df.symbol == symbol].cov_radius)

            self.cov_radii[symbol] = cov_radius
            return cov_radius
        
    def set_central_group(self, groupname):
        self.central_group_name = groupname

        # TODO make this variable
        groupnames = pd.read_csv("./data/central_groups.csv", header=0)

        df = groupnames[groupnames.name == groupname]

        self.central_group_atoms = {}

        for _, row in df.iterrows():
            self.central_group_atoms[row["atom"]] = int(row["amount"])

            if row.center == True and row.amount == 1:
                self.center_atom = row.atom
                self.center_ring = False
            elif row.center == True:
                self.center_atom = False
                self.center_ring = row.atom
        