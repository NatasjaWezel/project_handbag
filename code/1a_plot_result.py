from helpers.plot_functions import plot_fragments
from classes.Fragment import Fragment
from classes.Atom import Atom
from classes.Settings import Settings

from helpers.helpers import read_results_alignment

import pandas as pd

import sys
from tqdm import tqdm 


def main():
    if len(sys.argv) != 3:
        print("Usage: python run.py <inputfilename> <fragments_to_plot>")
        sys.exit(1)
    
    settings = Settings(sys.argv[1])
    df = read_results_alignment(settings.get_aligned_csv_filename())
    # df = df[df.in_central_group]

    amount = int(sys.argv[2])
   
    plot_fragments(df, amount)


if __name__ == "__main__":
    main()