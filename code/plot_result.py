from helpers.plot_functions import plot_fragments
from classes.Fragment import Fragment
from classes.Atom import Atom

from helpers.helpers import read_results_alignment

import pandas as pd

import sys
from tqdm import tqdm 


def main():
    if len(sys.argv) != 3:
        print("Usage: python run.py <inputfilename> <fragments_to_plot>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]
    amount = int(sys.argv[2])

    df = read_results_alignment(inputfilename)

    fragments = []

    unique_entrys = list(df.id.unique())[:amount]
    
    for unique_entry in tqdm(unique_entrys):
        atoms = df[df.id == unique_entry]
        
        fragment = Fragment(atoms.id.unique(), "")

        for _, row in atoms.iterrows():
            atom = Atom(row.atom_label, [row.atom_x, row.atom_y, row.atom_z])
            fragment.add_atom(atom)
        
        fragments.append(fragment)
    
    plot_fragments(fragments)


if __name__ == "__main__":
    main()