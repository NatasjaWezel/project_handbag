from classes.Atom import Atom
from classes.Fragment import Fragment

import pandas as pd
import math

import pickle

import time
from tqdm import tqdm

import csv

import copy

import numpy as np

def calculate_center(fragment_df):
    # TODO: this only works if there are no other C's on the ring
    atom_df = fragment_df[fragment_df.atom_symbol == 'C']

    coordinates = [atom_df.atom_x.mean(), atom_df.atom_y.mean(), atom_df.atom_z.mean()]

    assert (len(atom_df) == 6), "More then 6 C's in fragment, can't calculate center"

    return coordinates


def make_avg_fragment_if_not_exists(settings, df):
    try:
        fragment = pd.read_csv(settings.get_avg_fragment_filename())

        return fragment
    except FileNotFoundError:
        fragment = average_fragment(df)

        vdw_radii = [settings.get_vdw_radius(row.atom_symbol) for _, row in fragment.iterrows()]

        fragment['vdw_radius'] = vdw_radii
            
        fragment.to_csv(settings.get_avg_fragment_filename())

        return fragment


def average_fragment(df):
    """ Returns a fragment containing the average points of the central groups. """ 
    
    central_group_df = df[df.in_central_group]

    central_group_df = central_group_df.drop(columns=['entry_id', 'id', 'in_central_group'])
    avg_fragment_df = central_group_df.groupby('atom_label').agg({'atom_symbol': 'first', 'atom_x': 'mean', 'atom_y': 'mean', 'atom_z': 'mean'})


    return avg_fragment_df

def calculate_longest_vdw_radius_contact(fragment_df, settings):
    longest_distance = 0
    atom_a, atom_b = None, None

    # TODO: if there's an R, take the biggest vdw radius
    for _, atom1 in fragment_df.iterrows():
        for _, atom2 in fragment_df.iterrows():
            if not atom1.in_central_group and not atom2.in_central_group:
                distance = math.sqrt((atom1.atom_x - atom2.atom_x)**2 + (atom1.atom_y - atom2.atom_y)**2 + (atom1.atom_z - atom2.atom_z)**2)

                if distance > longest_distance:
                    longest_distance = distance
                    atom_a, atom_b = atom1, atom2

    longest_vdw_distance = (longest_distance + settings.get_vdw_radius(atom_a.atom_symbol) + settings.get_vdw_radius(atom_b.atom_symbol))/2

    return longest_vdw_distance
