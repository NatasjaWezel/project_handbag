from helpers import plot_fragments
import pandas as pd

from Atom import Atom
from Molecule import Fragment

def main():
    labels, fragments = read_csv(filename="results/saved_data1599472981.19713.csv")
    plot_fragments(fragments, labels=labels)

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