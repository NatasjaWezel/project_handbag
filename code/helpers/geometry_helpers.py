from classes.Atom import Atom
from classes.Fragment import Fragment

import pandas as pd
import math

import pickle

import time
from tqdm import tqdm

import csv

import copy

def calculate_center(fragment_df, atoms):
    frames = []

    for atom in atoms:
        atom_df = fragment_df[fragment_df.atom_symbol == atom]
        frames.append(atom_df)

    result = pd.concat(frames)
    coordinates = [result.atom_x.mean(), result.atom_y.mean(), result.atom_z.mean()]

    return coordinates


def make_avg_fragment_if_not_exists(settings, df):
    try:
        openpicklefile = open(settings.get_avg_fragment_filename(), 'rb')
        fragment = pickle.load(openpicklefile)
        openpicklefile.close()

        return fragment
    except FileNotFoundError:
        fragment = average_fragment(settings, df)
        
        openpicklefile = open(settings.get_avg_fragment_filename(), 'wb')
        pickle.dump(fragment, openpicklefile)
        openpicklefile.close()

        return fragment


def average_fragment(settings, df):
    """ Returns a fragment containing the average points of the central groups. """ 
    
    central_group_df = df[df.in_central_group]

    columns, atoms, count_dict = settings.get_avg_fragment_helpers()

    columns.append('time')

    # count how many atoms in one fragment
    new_df = pd.DataFrame(columns=columns, index=central_group_df.id.unique())

    # put coordinates in new df
    new_df = fill_coordinates(central_group_df, new_df, count_dict)

    # save it for stats
    new_df = new_df.apply(pd.to_numeric, downcast='float', errors='coerce')
    new_df.to_hdf(settings.get_avg_fragment_hdf_filename(), 'key')

    fragment = make_fragment(atoms, new_df)

    return fragment   


def fill_coordinates(central_group_df, new_df, count_dict):
    labels = central_group_df.id.unique()

    # put first fragment in there
    label = labels[0]
    single_fragment_df = central_group_df[central_group_df.id == label]

    closest = {}
    for _, row in single_fragment_df.iterrows():
        number = count_dict[row.atom_symbol]

        new_df["time"] = time.time()
        new_df = add_coordinates_to_df(new_df, label, row.atom_symbol + str(number), row)

        closest[row.atom_symbol + str(number)] = [row.atom_x, row.atom_y, row.atom_z]
        count_dict[row.atom_symbol] += 1

    print("Calculating average fragment: ")
    for label in tqdm(["BAPTIV1"]):
        single_fragment_df = central_group_df[central_group_df.id == label]
        closest_Hs = [i for i in closest if "H" in i]

        for _, row in single_fragment_df.iterrows():
            closest_atom, closest_Hs = get_closest_atom(closest_Hs, closest, row, new_df)
                
            new_df = add_coordinates_to_df(df=new_df, label=label, atom_symbol=closest_atom, row=row)
            
    dropped_na = new_df.dropna()
    print("Dropping", len(new_df) - len(dropped_na), "fragments for avg due to NaN values")
    return dropped_na

def get_closest_atom(closest_Hs, closest, row, new_df):
    print(row.atom_label, end="")
    new_df['time'] = time.time()
            
    distance = math.inf
    closest_atom = None

    # find which other atom in the central group it's closest too
    for key, value in closest.items():
        d = math.sqrt((row.atom_x - value[0])**2 + (row.atom_y - value[1])**2 + (row.atom_z - value[2])**2)

        if d < distance:
            distance = d
            closest_atom = key


    print(closest_Hs, end="")
    # TODO: what if there are two ch3 groups, can we make sure they get on either side
    # make sure that for a CH3 group for example all Hs get a "closest" H
    if "H" in closest_atom and not closest_atom in closest_Hs:
        closest_atom = closest_Hs.pop()
    elif "H" in closest_atom:
        closest_Hs.remove(closest_atom)
    
    print(closest_Hs, end="")
    print(closest_atom)
    return closest_atom, closest_Hs

def make_fragment(atoms, new_df):
    fragment = Fragment(from_entry="allentries", fragment_id=1)

    for atom_label in atoms:
        coordinates = [new_df[atom_label + "x"].mean(), new_df[atom_label + "y"].mean(), new_df[atom_label + "z"].mean()]

        atom = Atom(atom_label, coordinates)

        fragment.add_atom(atom) 
    
    return fragment

import numpy as np

def add_coordinates_to_df(df, label, atom_symbol, row):
    try:
        # this gives a keyerror if this row/column does not exist
        df.loc[df.index == label, atom_symbol + "x"]

        df.loc[df.index == label, atom_symbol + "x"] = row.atom_x
        df.loc[df.index == label, atom_symbol + "y"] = row.atom_y
        df.loc[df.index == label, atom_symbol + "z"] = row.atom_z
    except KeyError as e:
        print("Amount of nans:", df.loc[df.index == label].isnull().sum().sum())
        print(df.loc[df.index == label])

        assert int(df.loc[df.index == label].isnull().sum().sum()) == 3, "KeyError" + str(e) + "while calculating avg fragment: R is not the last atom, something else went wrong"

        df.loc[df.index == label, "R1x"] = row.atom_x
        df.loc[df.index == label, "R1y"] = row.atom_y
        df.loc[df.index == label, "R1z"] = row.atom_z

    return df
