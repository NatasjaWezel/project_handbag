import os

import pandas as pd
import math


class Settings():
    def __init__(self, WORKDIR, coordinate_file, labelfile):

        self.coordinate_data = coordinate_file
        self.label_data = labelfile

        names = coordinate_file.rsplit('\\')[-1].rsplit('.', 1)[0].rsplit('_aligned', 1)[0]

        self.central_group_name = names.split("_")[0]
        self.contact_group_name = names.split("_")[1]

        # setup results files
        self.output_folder_central_group = WORKDIR + "\\results\\" + names.split("_")[0] + "\\"
        output_folder_specific = self.output_folder_central_group + names + "\\"

        if not os.path.exists(self.output_folder_central_group):
            os.mkdir(self.output_folder_central_group)

        if not os.path.exists(output_folder_specific):
            os.mkdir(output_folder_specific)

        self.outputfile_prefix = output_folder_specific + names.rsplit('_', 1)[0]

    def get_directionality_results_filename(self):
        return self.output_folder_central_group + self.central_group_name + "_directionality_results.csv"

    def get_aligned_csv_filename(self):
        return self.outputfile_prefix + "_kabsch_aligned.csv"

    def get_coordinate_df_filename(self):
        return self.outputfile_prefix + "_coordinates_contact.hdf"

    def get_coordinate_df_key(self):
        return self.to_count_contact

    def get_density_df_filename(self):
        density_df_filename = self.outputfile_prefix + "_density.hdf"
        return density_df_filename

    def get_density_df_key(self):
        return self.to_count_contact + str(self.resolution).rstrip("0").replace(".", "")

    def set_resolution(self, resolution):
        self.resolution = resolution

    def get_density_plotname(self):
        density_plotname = self.outputfile_prefix + "_" + str(self.resolution) + "_density.svg"
        return density_plotname

    def set_atom_to_count(self, atom_str):
        self.to_count_contact = atom_str

    def get_avg_frag_filename(self):
        avg_fragment_filename = self.outputfile_prefix + "_avg_fragment.csv"
        return avg_fragment_filename


class Radii():
    def __init__(self, RADII_CSV):
        self.radii_filename = RADII_CSV
        self.vdw_radii = {}

    def get_vdw_radius(self, symbol):
        if symbol not in self.vdw_radii.keys():
            radii_df = pd.read_csv(self.radii_filename)

            vdw_radius = float(radii_df[radii_df.symbol == symbol].vdw_radius)

            self.vdw_radii[symbol] = vdw_radius

        return self.vdw_radii[symbol]

    def get_vdw_distance_contact(self, df, settings):
        if settings.to_count_contact == "centroid":
            return self.calculate_longest_vdw_radius_contact(df, settings)

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


class AlignmentSettings(Settings):
    def __init__(self, WORKDIR, coordinate_file, labelfile):
        Settings.__init__(self, WORKDIR, coordinate_file, labelfile)

        self.alignment = {}

    def set_central_group_csv(self, filename):
        self.central_group_csv = filename

    def get_central_group_csv(self):
        return self.central_group_csv

    def get_aligned_csv_filenames(self):
        aligned_csv = self.outputfile_prefix + "_aligned.csv"
        structures_csv = self.outputfile_prefix + "_structures.csv"
        return aligned_csv, structures_csv

    def prepare_alignment(self):
        self.read_coord_file()
        self.make_alignment_dict()

    def read_coord_file(self):
        """ Reads the first 100 lines of a csv file to count the atoms per fragment and per central
            group and which atom will get what label from the parameter file. """

        translation_dict = self.atoms_to_labels()

        with open(self.coordinate_data) as inputfile:
            lines = [next(inputfile) for x in range(100)]

        self.label_list = []

        # count number of atoms in a fragment
        self.no_atoms = 0

        # loop through the atoms until next fragment is found
        for line in lines[1:]:
            if "FRAG" in line:
                break

            # append to the label list which label the atom has, if any, else '-'
            self.label_list.append(translation_dict.get(line.split(' ')[0], '-'))
            self.no_atoms += 1

    def atoms_to_labels(self):
        # read the first two lines to figure out the label and atom ordeer
        with open(self.label_data) as inputfile:
            lines = [next(inputfile) for x in range(2)]

        # split it into two lists
        labels, atoms = lines[0].split(','), lines[1].split(',')

        # build a dictionary with atoms as keys and labels as their values
        translation_dict = {}
        for label, atom in zip(labels, atoms):
            if "LAB" in label:
                translation_dict[atom.strip()] = label.strip()

        return translation_dict

    def make_alignment_dict(self):
        """ """

        # get alignment info from file
        df = pd.read_csv(self.central_group_csv)
        df = df[df.name == self.central_group_name]

        # set up alignment dictionary
        self.alignment = {}
        self.alignment['center'] = df.center_label.max()
        self.alignment['yaxis'] = df.y_axis_label.max()
        self.alignment['xyplane'] = df.xy_plane_label.max()

        self.alignment['bin'] = df.bin.max()

        if self.alignment['bin'] != '-':
            self.alignment['bin'] = self.alignment['bin'].split('-')

        self.alignment['R'] = df.R.max() + '-' + df.not_R.max()

        # rename label for R group to 'Ri'
        if self.alignment['R'] != '-':
            R_atoms = self.alignment['R'].split('-')

            for i, Rlabel in enumerate(R_atoms):
                new_Rlabel = 'R' + str(i+1)
                self.label_list[self.label_list.index(Rlabel)] = new_Rlabel

                self.alignment[Rlabel] = Rlabel

                for key, value in self.alignment.items():
                    if value == Rlabel:
                        self.alignment[key] = new_Rlabel

    def get_index_alignment_atom(self, position):
        return self.label_list.index(self.alignment[position])
