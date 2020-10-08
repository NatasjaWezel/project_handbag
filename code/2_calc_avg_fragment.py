from helpers.geometry_helpers import make_avg_fragment_if_not_exists
from helpers.helpers import read_results_alignment
from classes.Settings import Settings

import pandas as pd
import numpy as np

import sys

def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_density.py <path/to/inputfile>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    settings = Settings(inputfilename)
    settings.set_central_group()

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())
    
    make_avg_fragment_if_not_exists(settings, aligned_fragments_df)


if __name__ == "__main__":
    main()