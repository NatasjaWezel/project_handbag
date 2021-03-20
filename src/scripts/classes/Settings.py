import os

import pandas as pd


class Settings():
    def __init__(self, WORKDIR, coordinate_file):
        self.WORKDIR = WORKDIR

        self.coordinate_file = coordinate_file

        name = coordinate_file.rsplit('\\')[-1].rsplit('.', 1)[0]
        print(name)

        self.central_name = name.split("_")[0]
        self.contact_name = name.split("_")[1]

        self.custom_alignment_file = False
        self.custom_structure_file = False
        self.threshold = 0.1
        self.resolution = False

        if 'data' in coordinate_file:
            # setup results files
            if not os.path.exists(WORKDIR +"\\results"):
                os.mkdir(WORKDIR + "\\results")
                os.mkdir(WORKDIR + "\\results\\pairs")
                
            self.output_folder_central_group = WORKDIR + "\\results\\pairs\\" + name.split("_")[0] + "\\"
            self.output_folder_specific = self.output_folder_central_group + name + "\\"

            if not os.path.exists(self.output_folder_central_group):
                os.mkdir(self.output_folder_central_group)
            self.outputfile_prefix = self.output_folder_specific + self.central_name + '_' + self.contact_name
        else:
            # setup results files
            self.output_folder_specific = coordinate_file.rsplit('\\', 1)[0] + '\\'
            self.outputfile_prefix = self.output_folder_specific + name

        if not os.path.exists(self.output_folder_specific):
            os.mkdir(self.output_folder_specific)

    def get_central_groups_csv_filename(self):
        return self.WORKDIR + "\\src\\data\\central_groups.csv"

    def get_methyl_csv_filename(self):
        return self.WORKDIR + "\\src\\data\\methylmodel.csv"

    def get_finger_print_filename(self):
        return self.WORKDIR + "\\src\\data\\fingerprints.csv"

    def get_radii_csv_name(self):
        return self.WORKDIR + "\\src\\data\\radii.csv"

    def get_directionality_results_filename(self):
        return self.output_folder_central_group + self.central_name + "_directionality_results.csv"

    def set_custom_alignment_filename(self, name):
        self.custom_alignment_file = self.output_folder_specific + name
        return self.custom_alignment_file

    def set_custom_structures_filename(self, name):
        self.custom_structure_file = self.output_folder_specific + name
        return self.custom_structure_file

    def get_aligned_csv_filename(self):
        if self.custom_alignment_file:
            return self.custom_alignment_file

        return self.outputfile_prefix + "_aligned.csv"

    def get_structure_csv_filename(self):
        if self.custom_structure_file:
            return self.custom_structure_file

        return self.outputfile_prefix + "_structures.csv"

    def get_coordinate_df_filename(self):
        return self.outputfile_prefix + "_coordinates_contact.hdf"

    def get_coordinate_df_key(self):
        return self.contact_rp

    def get_density_df_filename(self):
        density_df_filename = self.outputfile_prefix + "_density.hdf"
        return density_df_filename

    def get_density_df_key(self):
        return self.contact_rp + str(self.resolution).rstrip("0").replace(".", "")

    def set_resolution(self, resolution):
        self.resolution = resolution

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_density_plotname(self):
        density_plotname = self.outputfile_prefix + "_" + str(self.resolution) + "_density.svg"
        return density_plotname

    def set_contact_reference_point(self, atom_str):
        self.contact_rp = atom_str

    def get_avg_frag_filename(self):
        avg_fragment_filename = self.outputfile_prefix + "_avg_fragment.csv"
        return avg_fragment_filename


class AlignmentSettings(Settings):
    def __init__(self, WORKDIR, coordinate_file):
        Settings.__init__(self, WORKDIR, coordinate_file)

        self.label_data = coordinate_file.rsplit('.', 1)[0] + '.csv'
        self.alignment = {}

    def set_label_file(self, filename):
        self.label_data = filename

    def set_no_fragments(self, no_fragments):
        self.no_fragments = no_fragments

    def get_no_fragments(self):
        return self.no_fragments

    def get_aligned_csv_filenames(self):
        aligned_csv = self.get_aligned_csv_filename()
        structures_csv = self.get_structure_csv_filename()
        return aligned_csv, structures_csv

    def prepare_alignment(self):
        self.read_coord_file()
        self.make_alignment_dict()

    def update_coordinate_filename(self):
        # TODO: fix dit met multiple files and original file too big enzo
        if self.central_name == "RC6H5" and (self.contact_name == "ArCH" or
                                             self.contact_name == "RC6H5"):

            self.coordinate_file = self.coordinate_file.rsplit('.', 1)[0] + "_1."\
                                   + self.coordinate_file.rsplit('.', 1)[1]

    def read_coord_file(self):
        """ Reads the first 100 lines of a csv file to count the atoms per fragment and per central
            group and which atom will get what label from the parameter file. """

        translation_dict = self.atoms_to_labels()

        with open(self.coordinate_file) as inputfile:
            lines = [next(inputfile) for x in range(100)]

        self.label_list = []

        # count number of atoms in a fragment
        self.no_atoms = 0
        self.no_atoms_central = 0

        # loop through the atoms until next fragment is found
        for line in lines[1:]:
            if "FRAG" in line:
                break

            # append to the label list which label the atom has, if any, else '-'
            self.label_list.append(translation_dict.get(line.split(' ')[0], '-'))
            self.no_atoms += 1

            if translation_dict.get(line.split(' ')[0], '-') != '-':
                self.no_atoms_central += 1

    def atoms_to_labels(self):
        # read the first two lines to figure out the label and atom ordeer
        with open(self.label_data, 'r') as inputfile:
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
        df = pd.read_csv(self.get_central_groups_csv_filename())
        df = df[df.name == self.central_name]

        # set up alignment dictionary
        self.alignment = {}
        self.alignment['center'] = df.center_label.max()
        self.alignment['yaxis'] = df.y_axis_label.max()
        self.alignment['xyplane'] = df.xy_plane_label.max()

        self.alignment['bin'] = df.bin.max()

        if self.alignment['bin'] != '-':
            self.alignment['bin'] = self.alignment['bin'].split('-')

        self.alignment['R'] = df.R.max()
        self.alignment['r'] = df.not_R.max()

        # rename label for R group to 'Ri'
        if self.alignment['R'] != '-' and self.alignment['R'] != '':
            R_atoms = self.alignment['R'].split('-')

            for i, Rlabel in enumerate(R_atoms):
                new_Rlabel = Rlabel + '-R' + str(i + 1)
                self.label_list[self.label_list.index(Rlabel)] = new_Rlabel

                self.alignment[Rlabel] = Rlabel

                for key, value in self.alignment.items():
                    if value == Rlabel:
                        self.alignment[key] = new_Rlabel

        if self.alignment['r'] != '-' and self.alignment['r'] != '':
            R_atoms = self.alignment['r'].split('-')

            for i, Rlabel in enumerate(R_atoms):
                new_Rlabel = Rlabel + '-r' + str(i + 1)
                self.label_list[self.label_list.index(Rlabel)] = new_Rlabel

                self.alignment[Rlabel] = Rlabel

                for key, value in self.alignment.items():
                    if value == Rlabel:
                        self.alignment[key] = new_Rlabel

    def get_index_alignment_atom(self, position):
        return self.label_list.index(self.alignment[position])
