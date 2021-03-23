# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `Settings` is a class that contains almost all information about where to find datafiles, save result files, and find
# these result files again. It contains the names of the central and contact groups, and the specified threshold and
# resolution.
#
# It's child `AlignmentSettings` contains a bit more parameters specifically for superimposition.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
import sys

import pandas as pd


class Settings():
    """ Contains all information for where to find what file. Also contains central and contact names for name
        specification and parameters for threshold and resolution. """

    def __init__(self, WORKDIR, coordinate_file, central=False, contact=False):
        """ Init takes as input the Work directory so it can work from the main pipeline ánd from the scripts folder.
            It also takes the name of the coordinate file as input. When the name specifications are met, it can find
            all other necessary files itself. Else, the central and contact name have to be specified as well. """

        self.WORKDIR = WORKDIR
        self.coordinate_file = coordinate_file

        # setup results folders
        if not os.path.exists(WORKDIR + "\\results"):
            os.mkdir(WORKDIR + "\\results")
        if not os.path.exists(WORKDIR + "\\results\\pairs"):
            os.mkdir(WORKDIR + "\\results\\pairs")

        name = coordinate_file.rsplit('\\')[-1].rsplit('.', 1)[0]

        # check whether name specifications are met
        if central is False and contact is False:
            self.central_name = name.split("_")[0]
            self.contact_name = name.split("_")[1]

            self.set_result_directory(name)
        elif central is False or contact is False:
            print("Please specify both the central name and contact name.")
            sys.exit(1)
        else:
            self.central_name = central
            self.contact_name = contact

            # if the name is messy, give it a prefix containing central and group name
            self.set_result_directory(name, extra_prefix=central + "_" + contact + "_")

    def set_result_directory(self, name, extra_prefix=""):
        """ Set the results directory for all future generated files. """

        self.output_folder_central_group = self.WORKDIR + "\\results\\pairs\\" + self.central_name + "\\"
        self.output_folder_specific = self.output_folder_central_group + extra_prefix + name + "\\"

        self.outputfile_prefix = self.output_folder_specific + self.central_name + '_' + self.contact_name

        if not os.path.exists(self.output_folder_central_group):
            os.mkdir(self.output_folder_central_group)
        if not os.path.exists(self.output_folder_specific):
            os.mkdir(self.output_folder_specific)

    def set_resolution(self, resolution):
        self.resolution = round(resolution, 2)

    def set_threshold(self, threshold):
        self.threshold = round(threshold, 2)

    def get_aligned_csv_filename(self):
        return self.outputfile_prefix + "_aligned.csv"

    def get_structure_csv_filename(self):
        return self.outputfile_prefix + "_structures.csv"

    def set_contact_reference_point(self, atom_str):
        self.contact_rp = atom_str

    def get_central_groups_csv_filename(self):
        return self.WORKDIR + "\\src\\files\\central_groups.csv"

    def get_methyl_csv_filename(self):
        return self.WORKDIR + "\\src\\files\\methylmodel.csv"

    def get_fingerprint_filename(self):
        return self.WORKDIR + "\\src\\files\\fingerprints.csv"

    def get_radii_csv_name(self):
        return self.WORKDIR + "\\src\\files\\radii.csv"

    def get_directionality_results_filename(self):
        return self.output_folder_central_group + self.central_name + "_directionality_results.csv"

    def get_structure_csv_filename(self):
        return self.outputfile_prefix + "_structures.csv"

    def get_coordinate_df_filename(self):
        return self.outputfile_prefix + "_coordinates_contact.hdf"

    def get_coordinate_df_key(self):
        return self.contact_rp

    def get_density_df_filename(self):
        return self.outputfile_prefix + "_density.hdf"

    def get_density_df_key(self):
        return self.contact_rp + str(self.resolution).rstrip("0").replace(".", "")

    def get_density_plotname(self):
        return self.outputfile_prefix + "_" + str(self.resolution) + "_density.svg"

    def get_avg_frag_filename(self):
        avg_fragment_filename = self.outputfile_prefix + "_avg_fragment.csv"
        return avg_fragment_filename


class AlignmentSettings(Settings):
    """ Alignment Settings contains some extra parameters for superimposition, and inherits all functionality from its
        parent class Settings. All scripts that do superimposition have to use this class, otherwise only the settings
        class is necesarry. """

    def __init__(self, WORKDIR, coordinate_file, central=False, contact=False):
        Settings.__init__(self, WORKDIR, coordinate_file, central, contact)

        self.label_data = coordinate_file.rsplit('.', 1)[0] + '.csv'
        self.alignment = {}

        self.read_coord_file()
        self.make_alignment_dict()

    def set_label_file(self, filename):
        self.label_data = filename

    def set_no_fragments(self, no_fragments):
        self.no_fragments = no_fragments

    def get_no_fragments(self):
        return self.no_fragments

    def get_aligned_csv_filenames(self):
        return self.get_aligned_csv_filename(), self.get_structure_csv_filename()

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

        self.label_list = []

        # count number of atoms in a fragment
        self.no_atoms = 0
        self.no_atoms_central = 0

        with open(self.coordinate_file) as inputfile:

            # skip first row because that definitely contains "FRAG"
            next(inputfile)
            line = next(inputfile)

            # until the next fragment is found
            while "FRAG" not in line:

                # append to the label list which label the atom has, if any, else '-'
                self.label_list.append(translation_dict.get(line.split(' ')[0], '-'))
                self.no_atoms += 1

                if translation_dict.get(line.split(' ')[0], '-') != '-':
                    self.no_atoms_central += 1

                line = next(inputfile)

    def atoms_to_labels(self):
        """ Translate each atom its information from the label file to a label in the coordinate file. """

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
        """ Read the csv containing all the labels and make a dict containing new clean and recognizable labels, but
            also the old labels for backwards compatibility. """

        # get alignment info from file
        df = pd.read_csv(self.get_central_groups_csv_filename(), comment="#")
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
        self.alignment['r'] = df.treat_as_R.max()

        self.rename_labels('R', '-R')
        self.rename_labels('r', '-r')

    def rename_labels(self, label, new_label):
        """ Rename labels to recognizable as R group or treat as R group. """

        # rename label for LABx group to 'LABx-Ri'
        if self.alignment[label] != '-' and self.alignment[label] != '':
            R_atoms = self.alignment[label].split('-')

            # for each R label
            for i, Rlabel in enumerate(R_atoms):

                # make a new R label and update dict
                new_Rlabel = Rlabel + new_label + str(i + 1)
                self.label_list[self.label_list.index(Rlabel)] = new_Rlabel
                self.alignment[Rlabel] = Rlabel

                for key, value in self.alignment.items():
                    if value == Rlabel:
                        self.alignment[key] = new_Rlabel

    def get_index_alignment_atom(self, position):
        return self.label_list.index(self.alignment[position])
