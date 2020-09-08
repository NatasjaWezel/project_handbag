from helpers import plot_fragments
import pandas as pd

from Atom import Atom
from Molecule import Fragment

import pickle

from run import Save_molecules

def main():
    f = open('results/saved_data1599480077.3004494.pkl', 'rb')
    saved_molecules = pickle.load(f)
    f.close()

    fragments = []
    labels = []

    for molecule in saved_molecules.molecules:
        for fragment in molecule.fragments:
            fragments.append(fragment)
            labels.append(molecule.label + fragment.fragment_id)

    amount = 5
    plot_fragments(fragments[:amount], labels[:amount])


if __name__ == "__main__":
    main()