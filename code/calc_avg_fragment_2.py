import sys
import pandas as pd

from classes.Settings import Settings
from helpers.geometry_helpers import average_fragment, add_model_methyl


def main():

    if len(sys.argv) != 2:
        print("Usage: python 2_calc_avg_fragment.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    settings = Settings(inputfilename)

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())

    make_avg_fragment_if_not_exists(settings, aligned_fragments_df)


def read_results_alignment(inputfilename):
    df = pd.read_csv(inputfilename, header=None)
    df.columns = ["id", "entry_id", "atom_label", "atom_symbol", "in_central_group", "atom_x", "atom_y", "atom_z"]

    return df


def make_avg_fragment_if_not_exists(settings, df):
    try:
        fragment = pd.read_csv(settings.get_avg_fragment_filename())
    except FileNotFoundError:
        fragment = average_fragment(df, settings)

        print(settings.central_group_name)
        if settings.central_group_name == "RCOMe":
            fragment = add_model_methyl(fragment=fragment, settings=settings)

        fragment.to_csv(settings.get_avg_fragment_filename())

    return fragment


if __name__ == "__main__":
    main()
