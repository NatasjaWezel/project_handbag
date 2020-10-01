# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments. It then divides the 
# surrounding space into a number of bins, depending on which resolution is 
# set. It counts how many of the contact atoms/ centers of contact groups are
# are in each bin and normalizes that by the total amount of contact atoms or 
# groups. Then a plot is made that shows the density of the contacts in "4D". 
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.geometry_helpers import make_avg_fragment_if_not_exists
from helpers.helpers import read_results_alignment
from classes.Settings import Settings

import pandas as pd
import numpy as np

import sys



def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_density.py <path/to/inputfile> <central group>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    settings = Settings(inputfilename)
    settings.set_central_group(sys.argv[2])

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())
    
    make_avg_fragment_if_not_exists(settings, aligned_fragments_df)


if __name__ == "__main__":
    main()