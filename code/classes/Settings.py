import pandas as pd

from helpers.headers import RADII_CSV, CENTRAL_GROUPS_CSV, DATADIR, RESULTSDIR
import os

class Settings():
    def __init__(self, inputfilename):
        
        title = inputfilename.rsplit('\\')[-1].rsplit('.', 1)[0].rsplit('_aligned', 1)[0]

        self.central_group_name = title.split("_")[0]
        self.contact_group_name = title.split("_")[1]

        self.output_folder_central_group = RESULTSDIR + title.split("_")[0] + "\\"   
        output_folder_specific = self.output_folder_central_group + title + "\\"

        if not os.path.exists(self.output_folder_central_group):
            os.mkdir(self.output_folder_central_group)

        if not os.path.exists(output_folder_specific):
            os.mkdir(output_folder_specific)

        self.outputfile_prefix = output_folder_specific + title

        self.radii_filename = RADII_CSV

        self.vdw_radii = {}
        self.cov_radii = {}

    def get_directionality_results_filename(self):
        return self.output_folder_central_group + self.central_group_name + "_directionality_results.csv"

    def get_aligned_csv_filename(self):
        aligned_csv_name = self.outputfile_prefix + "_aligned.csv"
        return aligned_csv_name

    def get_current_tolerance(self):
        return self.tolerance

    def set_current_tolerance(self, tolerance):
        self.tolerance = tolerance
        
    def get_avg_fragment_filename(self):
        avg_fragment_filename = self.outputfile_prefix + "_avg_fragment.pkl"
        return avg_fragment_filename

    def get_avg_fragment_hdf_filename(self):
        avg_fragment_hdf_filename = self.outputfile_prefix + "_avg_fragment.hdf"
        return avg_fragment_hdf_filename

    def get_density_df_filename(self):
        density_df_filename = self.outputfile_prefix + "_density.hdf"
        return density_df_filename

    def set_resolution(self, resolution):
        self.resolution = resolution

    def get_density_plotname(self):
        density_plotname = self.outputfile_prefix + "_" + str(self.resolution) + "_density.pdf"
        return density_plotname

    def get_density_df_key(self):
        return "key" + str(self.resolution).rstrip("0").replace(".", "")

    def set_atom_to_count(self, atom_str):
        self.to_count_contact = atom_str

    def get_vdw_radius(self, symbol):
        
        if symbol in self.vdw_radii.keys():
            return self.vdw_radii[symbol]
        else:
            radii_df = pd.read_csv(self.radii_filename)

            vdw_radius = float(radii_df[radii_df.symbol == symbol].vdw_radius)

            self.vdw_radii[symbol] = vdw_radius
            return vdw_radius

    def get_avg_fragment_helpers(self):
        count_dict = {}

        for atom, amount in self.central_group_atoms.items():
            count_dict[atom] = amount

        return count_dict

    def get_cov_radius(self, symbol):
        
        if symbol in self.cov_radii.keys():
            return self.cov_radii[symbol]
        else:
            radii_df = pd.read_csv(self.radii_filename)

            cov_radius = float(radii_df[radii_df.symbol == symbol].cov_radius)

            self.cov_radii[symbol] = cov_radius
            return cov_radius
        
    def set_central_group(self):
        groupnames = pd.read_csv(CENTRAL_GROUPS_CSV, header=0)

        df = groupnames[groupnames.name == self.central_group_name]

        self.central_group_atoms = {}

        for _, row in df.iterrows():
            self.central_group_atoms[row["atom"]] = int(row["amount"])

            if row.center == True and row.amount == 1:
                self.center_atom = row.atom
                self.center_ring = False
            elif row.center == True:
                self.center_atom = False
                self.center_ring = row.atom

        self.amount_central_group_atoms = sum(self.central_group_atoms.values())
        