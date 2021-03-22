# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is the script that runs the entire pipeline of the program. As input, it needs a
# coordinate csv file and a parameter csv file that contains the labels of each atom of all
# fragments.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import argparse
import sys

import pandas as pd
import numpy as np

sys.path.append(".//scripts")

from classes.LoadArgsFromFile import LoadArgsFromFile
from classes.Settings import AlignmentSettings
from classes.Radii import Radii

from constants.paths import WORKDIR_MAIN
from constants.colors import COLORS

from helpers.plot_functions import plot_fragments

from align_kabsch import split_file_if_too_big, align_all_fragments
from plot_density import make_density_plot
from scripts.plot_avg_fragment_2a import plot_avg_fragment
from calc_avg_fragment import calc_avg_frag
from plot_contact_atoms import make_contact_rps_plot

from density_slider import make_density_slider_plot

from helpers.density_helpers import make_density_df, find_available_volume
from helpers.geometry_helpers import make_coordinate_df

from constants.constants import STANDARD_RES, STANDARD_THRESHOLD


def main():

    print("Welcome to the program!")
    parser = initiate_parser()
    args = check_args(parser)
    settings = make_settings_with_args(args)

    # Pipeline step 1: Align all fragments
    aligned_df = align_all_fragments(settings)

    # Pipeline step 2-3: Central group model
    radii = Radii(settings.get_radii_csv_name())
    central_model = calc_avg_frag(aligned_df, settings, radii)
    central_model.to_csv(settings.get_avg_frag_filename(), index=False)

    # Pipeline step 4: Distance contact atom/center to the model
    coordinate_df = make_coordinate_df(aligned_df, settings, central_model, radii)

    # Pipeline step 5: Density calculation
    density_df = make_density_df(settings, coordinate_df)
    density_df['datafrac_normalized'] = density_df[settings.contact_rp] / density_df[settings.contact_rp].sum()

    # Pipeline step 6: Volumes
    tolerance = 0.5
    contact_group_radius = radii.get_vdw_distance_contact(settings.contact_rp)
    Vavailable = find_available_volume(avg_fragment=central_model, extra=(tolerance + contact_group_radius))
    threshold_calc = density_df.datafrac_normalized.max() * settings.threshold
    in_cluster = density_df[density_df.datafrac_normalized >= threshold_calc]

    datafrac = in_cluster.datafrac_normalized.sum()
    Vcluster = len(in_cluster) * settings.resolution**3

    # Pipeline step 7: Directionality
    directionality = datafrac / Vcluster * (Vavailable/2)
    print(f"The directionality of {settings.central_name}--{settings.contact_name} ({settings.contact_rp}) is\
            {directionality}\n")

    # when done running, give option menu
    print_menu()
    possible_inputs = [1, 2, 3, 4, 5, 6, 7]
    option = ask_int_input("What do you want to plot?", possible_inputs)
    while not option == 7:
        perform_option(option, settings)
        print_menu()
        option = ask_int_input("What do you want to plot?", possible_inputs)

    print_epilog()


def perform_option(option, settings):
    # do something
    if option == 1:
        data = pd.read_csv(settings.get_aligned_csv_filename())
        max_frags = len(data.fragment_id.unique())
        possible_inputs = range(0, max_frags + 1)
        amount = ask_int_input("How many superimposed fragments would you like to plot?\n(Recommended < 100)\n",
                               possible_inputs)

        default = "Y"
        only_central = ask_bool_input("Do you want to plot the contact groups as well? [Y]\\N\n", default)
        print()

        if "y" == only_central.lower():
            data = data[data.label != "-"]
        plot_fragments(data, amount, COLORS)
    elif option == 2:
        plot_avg_fragment(settings)
    elif option == 3:
        pass
    elif option == 4:
        avg_fragment = pd.read_csv(settings.get_avg_frag_filename())
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
        make_density_plot(avg_fragment, density_df, settings)
    elif option == 5:
        df = pd.read_csv(settings.get_aligned_csv_filename())

        avg_fragment = pd.read_csv(settings.get_avg_frag_filename())

        radii = Radii(settings.get_radii_csv_name())
        coordinate_df = make_coordinate_df(df, settings, avg_fragment, radii)

        make_contact_rps_plot(avg_fragment, coordinate_df, settings)
    elif option == 6:
        default = "Y"
        confimation = ask_bool_input("The program still has to calculate non-standard resolutions." +
                                     "This may take some time. Do you want to continue? [Y]\\N\n", default)
        print()
        if (confimation.lower() == "y"):
            alligned_df = pd.read_csv(settings.get_aligned_csv_filename())
            central_model = pd.read_csv(settings.get_avg_frag_filename())
            radii = Radii(settings.get_radii_csv_name())
            coordinate_df = make_coordinate_df(alligned_df, settings, central_model, radii)
            for res in np.arange(0.2, 1.05, 0.05):
                settings.set_resolution(res)
                make_density_df(settings, coordinate_df)

            settings.set_resolution(STANDARD_RES)
            avg_fragment = pd.read_csv(settings.get_avg_frag_filename())
            make_density_slider_plot(avg_fragment, settings)


def ask_bool_input(message, default):
    option = input(message)
    valid = False
    while(not valid):
        try:
            if (option.lower() == "y" or option.lower() == "n"):
                valid = True
            elif (option.lower() == "yes" or option.lower() == "no"):
                valid = True
                option = option[0]
            elif option == "":
                option = default
                valid = True
            else:
                raise ValueError
        except ValueError:
            if (option == "quit") or (option == "exit"):
                print_epilog()
                exit()
            option = input("Please type Y or N. \n")

    return option


def ask_int_input(message, possible_inputs):
    option = input(message + "\n")

    valid = False
    while(not valid):
        try:
            val = int(option)
            if val in possible_inputs:
                valid = True
            else:
                raise ValueError
        except ValueError: 
            if (option == "quit") or (option == "exit"):
                print_epilog()
                exit()
            option = input("Type the integer belonging to your choice.\n")

    print()
    return int(option)


def print_epilog():
    print("Thank you for using our program!\nNatasja Wezel, Tiddo Mooibroek\nE-mail your suggestions, remarks and/or "
          + "questions to natasjawezel@gmail.com\nFind the source code at www.github.com/NatasjaWezel/MasterProject")


def print_menu():
    print("What plot do you want to see?")
    print("Plots")
    print("\t1 - aligned fragments")
    print("\t2 - model central group")
    print("\t3 - coordinates contact group")
    print("\t4 - density around central group")
    print("Interactive plots")
    print("\t5 - coordinates contact group")
    print("\t6 - density around central group")
    print("7 - to end the program")
    print()


def make_settings_with_args(args):

    if args.central is None and args.contact is None:
        settings = AlignmentSettings(WORKDIR=WORKDIR_MAIN, coordinate_file=args.input)
    if (args.central is not None and args.contact is None) or (args.central is None and args.contact is not None):
        print("You can only specify no contact/central group, or specify both..")
        sys.exit(1)
    if args.central is not None and args.contact is not None:
        settings = AlignmentSettings(WORKDIR=WORKDIR_MAIN, coordinate_file=args.input, central=args.central, contact=args.contact)
    if args.labels is not None:
        settings.set_label_file(args.labels)
    if args.vanderwaals is not None:
        pass
    if args.output is not None:
        pass

    settings.set_contact_reference_point(args.contact_rp.upper())
    settings.set_resolution(STANDARD_RES)
    settings.set_threshold(STANDARD_THRESHOLD)

    settings.prepare_alignment()

    split_file_if_too_big(settings.coordinate_file, settings.no_atoms)
    settings.update_coordinate_filename()

    print(f"Find your results in the output folder: {settings.output_folder_central_group}\n")

    return settings


def check_args(parser):
    args = parser.parse_args()

    if args.input is None and args.contact_rp is None:
        print('quantify.py: error: the following arguments are required: -i/--input, -crp/--contact_rp')
        sys.exit(1)
    elif args.contact_rp is None:
        print('quantify.py: error: the following arguments are required: -crp/--contact_rp')
        sys.exit(1)
    elif args.input is None:
        print('quantify.py: error: the following arguments are required: -i/--input, -crp/--contact_rp')
        sys.exit(1)

    return args


def initiate_parser():
    parser = argparse.ArgumentParser()

    # pop required group from parser to add one ourselves
    parser._action_groups.pop()
    # add required group first so it appears above optional
    required = parser.add_argument_group('required arguments')

    required.add_argument('-i', '--input', help='data input file from conquest')
    required.add_argument('-crp', '--contact_rp', help='atom or center of contact group of which we want the density')

    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('-f', '--argsfile', type=open, action=LoadArgsFromFile, help='load specified arguments from\
                          file')
    optional.add_argument('-l', '--labels', help='data input labels from conquest (default same as inputfilename but\
                          then .csv instead of .cor)')
    optional.add_argument('-ce', '--central', help='name of the central group (default first argument of inputfilename,\
                          splitted on underscores)')
    optional.add_argument('-co', '--contact', help='name of the contact group (default second argument of inputfilename,\
                          splitted on underscores)')
    optional.add_argument('-o', '--output', help='prefix of the outputfile (default CENTRAL_CONTACT)')
    optional.add_argument('-vdw', '--vanderwaals', help='used van der waals tolerance (default 0.5)')

    return parser


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_epilog()
