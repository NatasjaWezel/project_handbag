from classes.Atom import Atom
from classes.Fragment import Fragment

import pandas as pd
import math

import pickle

import time
from tqdm import tqdm

import csv

import copy

def calculate_center(fragment_df):
    #TODO: this only works if there are no other C's on the ring
    atom_df = fragment_df[fragment_df.atom_symbol == 'C']

    coordinates = [atom_df.atom_x.mean(), atom_df.atom_y.mean(), atom_df.atom_z.mean()]

    assert (len(atom_df) == 6), "More then 6 C's in fragment, can't calculate center"

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

    count_dict = settings.get_avg_fragment_helpers()

    atoms, columns = get_columns(central_group_df, count_dict)

    new_df = pd.DataFrame(columns=columns, index=central_group_df.id.unique())

    # put coordinates in new df
    new_df = fill_coordinates(central_group_df, new_df, atoms)

    # save it for stats
    new_df = new_df.apply(pd.to_numeric, downcast='float', errors='coerce')
    new_df.to_hdf(settings.get_avg_fragment_hdf_filename(), 'key')

    fragment = make_fragment(atoms, new_df)

    return fragment   

def get_columns(central_group_df, count_dict):
    label1 = central_group_df.id.unique()[0]
    single_fragment_df = central_group_df[central_group_df.id == label1]
    atom_order_df = list(single_fragment_df.atom_symbol)

    atoms = []
    amount_R = count_dict.get("R", 0)

    temp = {}
    for key in count_dict.keys():
        temp[key] = 1

    for atom in atom_order_df[:len(atom_order_df) - amount_R]:
        atoms.append(atom + str(temp[atom]))
        temp[atom] += 1
    
    for i in range(1, amount_R + 1):
        atoms.append("R" + str(i))

    cols = []
    for i in atoms:
        cols.extend([i + "x", i + "y", i + "z"])

    cols.append('time')

    return atoms, cols

def fill_coordinates(central_group_df, new_df, atoms):
    labels = central_group_df.id.unique()

    print("Calculating average fragment: ")
    for label in tqdm(labels):
        single_fragment_df = central_group_df[central_group_df.id == label]

        new_df.loc[new_df.index == label, "time"] = time.time()

        for i, row in enumerate(single_fragment_df.iterrows()):
            row = row[1]
            new_df = add_coordinates_to_df(df=new_df, label=label, row=row, atom=atoms[i])
            
    dropped_na = new_df.dropna()
    print("Dropping", len(new_df) - len(dropped_na), "fragments for avg due to NaN values")
    return dropped_na


def make_fragment(atoms, new_df):
    fragment = Fragment(from_entry="allentries", fragment_id=1)

    for atom_label in atoms:
        coordinates = [new_df[atom_label + "x"].mean(), new_df[atom_label + "y"].mean(), new_df[atom_label + "z"].mean()]

        atom = Atom(atom_label, coordinates)

        fragment.add_atom(atom) 
    
    print(fragment)
    return fragment

def add_coordinates_to_df(df, label, row, atom):

    df.loc[df.index == label, atom + "x"] = row.atom_x
    df.loc[df.index == label, atom + "y"] = row.atom_y
    df.loc[df.index == label, atom + "z"] = row.atom_z

    return df

def calculate_longest_vdw_radius_contact(fragment_df, settings):
    longest_distance = 0
    atom_a, atom_b = None, None

    for _, atom1 in fragment_df.iterrows():
        for _, atom2 in fragment_df.iterrows():
            if not atom1.in_central_group and not atom2.in_central_group:
                distance = math.sqrt((atom1.atom_x - atom2.atom_x)**2 + (atom1.atom_y - atom2.atom_y)**2 + (atom1.atom_z - atom2.atom_z)**2)

                if distance > longest_distance:
                    longest_distance = distance
                    atom_a, atom_b = atom1, atom2

    longest_vdw_distance = (longest_distance + settings.get_vdw_radius(atom_a.atom_symbol) + settings.get_vdw_radius(atom_b.atom_symbol))/2

    return longest_vdw_distance
