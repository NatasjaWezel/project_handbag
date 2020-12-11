import sys
import pandas as pd

from classes.Settings import Settings
from helpers.geometry_helpers import average_fragment, add_model_methyl


def main():

    if len(sys.argv) != 2:
        print("Usage: python calc_avg_fragment_2.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    avg_frag_settings = Settings(inputfilename)

    df = pd.read_csv(avg_frag_settings.get_kabsch_aligned_csv_filename(), header=0)

    avg_frag_file = avg_frag_settings.get_avg_frag_filename()

    fragment = average_fragment(df, avg_frag_settings)

    if avg_frag_settings.central_group_name == "RCOMe":
        fragment = add_model_methyl(fragment=fragment, settings=avg_frag_settings)

    fragment.to_csv(avg_frag_file, index=False)

    return fragment


if __name__ == "__main__":
    main()
