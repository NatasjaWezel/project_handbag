from classes.Atom import Atom
from classes.Fragment import Fragment

import pandas as pd
import math

import pickle

from tqdm import tqdm

def calculate_center(fragment_df, atoms):
    frames = []

    for atom in atoms:
        atom_df = fragment_df[fragment_df.atom_symbol == atom]
        frames.append(atom_df)

    result = pd.concat(frames)
    coordinates = [result.atom_x.mean(), result.atom_y.mean(), result.atom_z.mean()]

    return coordinates

def average_fragment(settings, df):
    """ Returns a fragment containing the bonds and average points of the interacting central groups. 
        # TODO: its NO3 specific right now """ 

    try:
        openpicklefile = open(settings.get_avg_fragment_filename(), 'rb')
        fragment = pickle.load(openpicklefile)
        openpicklefile.close()
    except FileNotFoundError:
        central_group_df = df[df.in_central_group]

        columns = []
        atoms = []

        count_dict = {}

        for atom, amount in settings.central_group_atoms.items():
            count_dict[atom] = 1

            for i in range(1, amount + 1):
                atoms.append(atom + str(i))
                columns.extend([atom + str(i) + "x", atom + str(i) + "y", atom + str(i) + "z"])

        # count how many atoms in one fragment
        new_df = pd.DataFrame(columns=columns, index=central_group_df.id.unique())

        # put first fragment in there
        label = central_group_df.id.unique()[0]
        single_fragment_df = central_group_df[central_group_df.id == label]

        closest = {}
        for _, row in single_fragment_df.iterrows():
            number = count_dict[row.atom_symbol]

            new_df = add_coordinates_to_df(new_df, label, row.atom_symbol + str(number), row)

            closest[row.atom_symbol + str(number)] = [row.atom_x, row.atom_y, row.atom_z]

            count_dict[row.atom_symbol] += 1

        # put rest of fragments in the new df
        new_df = fill_coordinates(central_group_df, new_df, closest)
        
        fragment = make_fragment(atoms, new_df)

        openpicklefile = open(settings.get_avg_fragment_filename(), 'wb')
        pickle.dump(fragment, openpicklefile)
        openpicklefile.close()
            
    return fragment   


def fill_coordinates(central_group_df, new_df, closest):
    # TODO: time this and check std's to see what's worth and what's not
    labels = central_group_df.id.unique()[1:]
    
    print("Calculating average fragment: ")
    for label in tqdm(labels):
        single_fragment_df = central_group_df[central_group_df.id == label]

        for _, row in single_fragment_df.iterrows():
            distance = math.inf
            closest_atom = None

            # find which O it's closest to
            for key, value in closest.items():
                d = math.sqrt((row.atom_x - value[0])**2 + (row.atom_y - value[1])**2 + (row.atom_z - value[2])**2)

                if d < distance:
                    distance = d
                    closest_atom = key

            new_df = add_coordinates_to_df(df=new_df, label=label, atom_symbol=closest_atom, row=row)

    return new_df

def make_fragment(atoms, new_df):
    fragment = Fragment(from_entry="allentries", fragment_id=1)

    for atom_label in atoms:
        coordinates = [new_df[atom_label + "x"].mean(), new_df[atom_label + "y"].mean(), new_df[atom_label + "z"].mean()]

        atom = Atom(atom_label, coordinates)

        fragment.add_atom(atom) 
    
    return fragment


def add_coordinates_to_df(df, label, atom_symbol, row):
    if label == "FOYMIS1":
        print(df.loc[df.index == label], atom_symbol)
    df.loc[df.index == label, atom_symbol + "x"] = row.atom_x
    df.loc[df.index == label, atom_symbol + "y"] = row.atom_y
    df.loc[df.index == label, atom_symbol + "z"] = row.atom_z

    return df
