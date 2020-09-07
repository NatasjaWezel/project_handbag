from helpers import plot_fragments
import pandas as pd

from Atom import Atom
from Molecule import Fragment

import pickle

from run import Save_molecules

def main():
    # labels, fragments = read_csv(filename="results/saved_data1599472981.19713.csv")
    # plot_fragments(fragments, labels=labels)

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

def read_csv(filename):
    df = pd.read_csv(filename, header=None)

    fragments = []
    labels = []

    df.columns = ['entry_id', 'fragment_id', 'atom_label', 'atom_x', 'atom_y', 'atom_z']

    entries = df.entry_id.unique()

    for entry in entries:
        smalldf = df[df.entry_id == entry]

        fragment_ids = smalldf.fragment_id.unique()

        for fragment_id in fragment_ids:
            fragment_df = smalldf[smalldf.fragment_id == fragment_id]

            fragment = read_fragment_from_csv(fragment_df)
            fragments.append(fragment)
            labels.append("Entry: " + entry + " fragment id: " + str(fragment_id))

    return labels, fragments

def read_fragment_from_csv(df):
    fragment = Fragment(int(df.fragment_id.unique()))

    for row in df.iterrows():
        items = row[1]
        # print(items.atom_label)
        atom = Atom(items.atom_label, items.atom_label, items.atom_x, items.atom_y, items.atom_z)
        fragment.atoms[items.atom_label] = atom

    return fragment


if __name__ == "__main__":
    main()