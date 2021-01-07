import sys
import pandas as pd

from classes.Settings import Settings, Radii
from helpers.geometry_helpers import average_fragment, add_model_methyl

from constants.paths import WORKDIR, RADII_CSV


def main():

    if len(sys.argv) != 2:
        print("Usage: python calc_avg_fragment_2.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    avg_frag_settings = Settings(WORKDIR, inputfilename)

    df = pd.read_csv(avg_frag_settings.get_aligned_csv_filename(), header=0)

    # make radii object to get vdw radii
    radii = Radii(RADII_CSV)
    fragment = average_fragment(df, avg_frag_settings, radii)

    if avg_frag_settings.central_group_name == "RCOMe":
        fragment = add_model_methyl(fragment=fragment, settings=avg_frag_settings)

    # get name and save average fragment
    avg_frag_file = avg_frag_settings.get_avg_frag_filename()
    fragment.to_csv(avg_frag_file, index=False)

    return fragment


if __name__ == "__main__":
    main()
