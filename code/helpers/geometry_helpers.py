from helpers.Atom import Atom
from helpers.Fragment import Fragment

import pandas as pd
import math

import pickle

def calculate_center(fragment_df, atoms):
    frames = []

    for atom in atoms:
        atom_df = fragment_df[fragment_df.atom_symbol == atom]
        frames.append(atom_df)

    result = pd.concat(frames)
    coordinates = [result.atom_x.mean(), result.atom_y.mean(), result.atom_z.mean()]

    return coordinates

def average_fragment(avg_fragment_name, df):
    """ Returns a fragment containing the bonds and average points of the interacting central groups. 
        # TODO: its NO3 specific right now """ 

    try:
        openpicklefile = open(avg_fragment_name, 'rb')
        fragment = pickle.load(openpicklefile)
        openpicklefile.close()
    except FileNotFoundError:
        df["unique_f_label"] = df["entry_id"] + df["fragment_id"].astype(str)

        central_group_df = df[df.fragment_or_contact == "c"]

        ideal_atoms = ["N", "O1", "O2", "O3"]

        columns = []
        for atom in ideal_atoms:
            columns.extend([atom + "x", atom + "y", atom + "z"])

        # count how many atoms in one fragment
        new_df = pd.DataFrame(columns=columns, index=central_group_df.unique_f_label.unique())

        # put first fragment in there
        label = central_group_df.unique_f_label.unique()[0]
        single_fragment_df = central_group_df[central_group_df.unique_f_label == label]

        counter = 1
        closest = {}

        for _, row in single_fragment_df.iterrows():
            if row.atom_symbol == "N":
                new_df = add_coordinates_to_df(new_df, label, row.atom_symbol, row)
            else:
                atom_symbol = "O" + str(counter)
                new_df = add_coordinates_to_df(new_df, label, atom_symbol, row)

                closest[atom_symbol] = [row.atom_x, row.atom_y, row.atom_z]
                counter += 1

        

        # TODO: time this and check std's to see what's worth and what's not
        labels = central_group_df.unique_f_label.unique()[1:]

        for label in labels:
            single_fragment_df = central_group_df[central_group_df.unique_f_label == label]

            for _, row in single_fragment_df.iterrows():
                if row.atom_symbol == "N":
                    new_df = add_coordinates_to_df(df=new_df, label=label, atom_symbol=row.atom_symbol, row=row)
                else:
                    distance = math.inf
                    closest_atom = None

                    # find which O it's closest to
                    for key, value in closest.items():
                        d = math.sqrt((row.atom_x - value[0])**2 + (row.atom_y - value[1])**2 + (row.atom_z - value[2])**2)

                        if d < distance:
                            distance = d
                            closest_atom = key

                    new_df = add_coordinates_to_df(df=new_df, label=label, atom_symbol=closest_atom, row=row)
        
        fragment = make_fragment(ideal_atoms, new_df)

        fragment.set_vdw_radii("data/vdw_radii.csv")

        openpicklefile = open(avg_fragment_name, 'wb')
        pickle.dump(fragment, openpicklefile)
        openpicklefile.close()
            
    return fragment   


def make_fragment(ideal_atoms, new_df):
    fragment = Fragment(from_entry="allentries", fragment_id=1)

    for atom_label in ideal_atoms:
        coordinates = [new_df[atom_label + "x"].mean(), new_df[atom_label + "y"].mean(), new_df[atom_label + "z"].mean()]

        atom = Atom(atom_label, coordinates)

        fragment.add_atom(atom) 
    
    return fragment


def add_coordinates_to_df(df, label, atom_symbol, row):
    df.loc[df.index == label, atom_symbol + "x"] = row.atom_x
    df.loc[df.index == label, atom_symbol + "y"] = row.atom_y
    df.loc[df.index == label, atom_symbol + "z"] = row.atom_z

    return df
