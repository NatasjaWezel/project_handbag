import sys

from classes.Settings import Settings
from helpers.geometry_helpers import make_avg_fragment_if_not_exists
from helpers.helpers import read_results_alignment


def main():

    if len(sys.argv) != 2:
        print("Usage: python 2_calc_avg_fragment.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    settings = Settings(inputfilename)

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())

    make_avg_fragment_if_not_exists(settings, aligned_fragments_df)


if __name__ == "__main__":
    main()
