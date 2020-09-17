from helpers.Atom import Atom
from helpers.Fragment import Fragment

def calculate_center(fragment_df, atoms):
    # TODO: only count "C" (or other atoms in atoms) for average
    atom_df = fragment_df
    x, y, z = atom_df.atom_x.mean(), atom_df.atom_y.mean(), atom_df.atom_z.mean()

    return x, y, z

def average_molecule(df):
    """ Returns a fragment containing the bonds and average points of the interacting molecule. """ 

    # TODO: make it an actual average
    # TODO: know which O is which O?
    central_group_df = df[df.fragment_or_contact == "c"]

    entry_id = central_group_df.entry_id.unique()[0]
    fragment_id = central_group_df.fragment_id.unique()[0]

    fragment = Fragment(entry_id, fragment_id)

    for _, row in central_group_df[(central_group_df.fragment_id == fragment_id) & (central_group_df.entry_id == entry_id)].iterrows():
        atom = Atom(row.atom_label, row.atom_x, row.atom_y, row.atom_z)
        fragment.add_atom(atom)

    return fragment    
