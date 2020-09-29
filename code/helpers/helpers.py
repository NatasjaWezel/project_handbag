from collections import defaultdict

import pandas as pd

def read_results_alignment(inputfilename):
    df = pd.read_csv(inputfilename, header=None)
    df.columns = ["id", "entry_id", "atom_label", "atom_symbol", "in_central_group", "atom_x", "atom_y", "atom_z"]

    return df

def check_if_label_exists(atom, fragment):
    """ Checks if label already exists in the fragment. Adds the laetter 'a' to it if it does. """
        
    if atom.label in fragment.atoms.keys():
        atom.label += 'a'
        atom = check_if_label_exists(atom, fragment)
    
    return atom
