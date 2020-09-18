from helpers.Atom import Atom
from helpers.Fragment import Fragment

import pandas as pd
import math

def calculate_center(fragment_df, atoms):
    # TODO: only count "C" (or other atoms in atoms) for average
    atom_df = fragment_df
    coordinates = [atom_df.atom_x.mean(), atom_df.atom_y.mean(), atom_df.atom_z.mean()]

    return coordinates

def average_fragment(df):
    """ Returns a fragment containing the bonds and average points of the interacting central groups. 
        # TODO: its NO3 specific right now """ 

    df["unique_f_label"] = df["entry_id"] + df["fragment_id"].astype(str)

    central_group_df = df[df.fragment_or_contact == "c"]

    ideal_atoms = ["N", "O1", "O2", "O3"]

    # count how many atoms in one fragment
    new_df = pd.DataFrame(columns=['Nx', 'Ny', 'Nz',
                                        'O1x', 'O1y', 'O1z', 
                                        'O2x', 'O2y', 'O2z', 
                                        'O3x', 'O3y', 'O3z'], index=central_group_df.unique_f_label.unique())

    # put first fragment in there
    label = central_group_df.unique_f_label.unique()[0]
    single_fragment_df = central_group_df[central_group_df.unique_f_label == label]

    O_counter = 1

    closest = {}

    # TODO: count the columnames here too
    for _, row in single_fragment_df.iterrows():
        if row.atom_element == "N":
            new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "Nx"] = row.atom_x
            new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "Ny"] = row.atom_y
            new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "Nz"] = row.atom_z
        else:
            new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "O" + str(O_counter)  + "x"] = row.atom_x
            new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "O" + str(O_counter)  + "y"] = row.atom_y
            new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "O" + str(O_counter)  + "z"] = row.atom_z
            closest["O" + str(O_counter)] = [row.atom_x, row.atom_y, row.atom_z]
            O_counter += 1

    # TODO: time this and check std's to see what's worth and what's not
    labels = central_group_df.unique_f_label.unique()[1:]

    for label in labels:
        single_fragment_df = central_group_df[central_group_df.unique_f_label == label]

        for _, row in single_fragment_df.iterrows():
            if row.atom_element == "N":
                new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "Nx"] = row.atom_x
                new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "Ny"] = row.atom_y
                new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], "Nz"] = row.atom_z
            else:
                distance = math.inf
                closest_atom = None

                # find which O it's closest to
                for key, value in closest.items():
                    d = math.sqrt((row.atom_x - value[0])**2 + (row.atom_y - value[1])**2 + (row.atom_z - value[2])**2)

                    if d < distance:
                        distance = d
                        closest_atom = key

                new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], closest_atom  + "x"] = row.atom_x
                new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], closest_atom  + "y"] = row.atom_y
                new_df.loc[new_df.index == single_fragment_df.unique_f_label.unique()[0], closest_atom  + "z"] = row.atom_z
    
    fragment = Fragment(from_entry="allentries", fragment_id=1)

    for atom_label in ideal_atoms:
        x = new_df[atom_label + "x"].mean()
        y = new_df[atom_label + "y"].mean()
        z = new_df[atom_label + "z"].mean()

        coordinates= [x,y,z]

        atom = Atom(atom_label, coordinates)

        fragment.add_atom(atom)
    
    print(fragment)
    
    return fragment    
