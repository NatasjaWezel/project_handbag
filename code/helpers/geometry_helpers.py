import math

import numpy as np
import pandas as pd
from tqdm import tqdm

import copy


def make_coordinate_df(df, settings, avg_fragment):
    try:
        coordinate_df = pd.read_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

        return coordinate_df
    except (KeyError, FileNotFoundError):
        first_fragment_df = df[df.id == df.id.unique()[0]]

        if settings.to_count_contact == "centroid":
            # plot centroids of all contact fragments
            coordinate_df = df.groupby("id").mean()
        elif len(first_fragment_df[first_fragment_df["atom_symbol"] == settings.to_count_contact]) == 1:
            # atom is unique, plot all of them
            coordinate_df = df[df.atom_symbol == settings.to_count_contact]
        else:
            # TODO: atom is not unique, find closest
            pass

        coordinate_df = distances_closest_vdw_central(coordinate_df, avg_fragment, settings)

        coordinate_df.to_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

        return coordinate_df


def distances_closest_vdw_central(coordinate_df, avg_fragment, settings):
    closest_distances = []
    closest_atoms_vdw = []

    points_avg_f = np.array([avg_fragment.atom_x, avg_fragment.atom_y, avg_fragment.atom_z]).T

    print("Searching for nearest atom from contact group...")
    print(len(coordinate_df.atom_x))
    for x, y, z in tqdm(zip(coordinate_df.atom_x, coordinate_df.atom_y, coordinate_df.atom_z)):

        p2 = np.array([x, y, z])

        dist = np.sqrt([np.sum((f - p2)**2, axis=0) for f in points_avg_f])

        min_dist_idx = dist.argmin()
        min_dist = dist[min_dist_idx]

        min_atom_vdw = avg_fragment.iloc[min_dist_idx]['vdw_radius']

        closest_distances.append(min_dist)
        closest_atoms_vdw.append(min_atom_vdw)

    coordinate_df["distance"] = closest_distances
    coordinate_df["vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df


def make_avg_fragment_if_not_exists(settings, df):
    try:
        fragment = pd.read_csv(settings.get_avg_fragment_filename())

        return fragment
    except FileNotFoundError:
        fragment = average_fragment(df, settings)

        if settings.central_group_name == "RCOMe":
            add_model_methyl(fragment)

        fragment.to_csv(settings.get_avg_fragment_filename())


def add_model_methyl(fragment, settings):
    # TODO: use labels

    print("Adding model CH3 group")

    # TODO: drop old H's not hardcoded
    fragment = fragment[(fragment.atom_label != "H5") & (fragment.atom_label != "H6") &
                        (fragment.atom_label != "H7")]

    a = np.array([float(fragment[fragment.index == "O3"].atom_x),
                  float(fragment[fragment.index == "O3"].atom_y),
                  float(fragment[fragment.index == "O3"].atom_z)])

    b = np.array([float(fragment[fragment.index == "C2"].atom_x),
                  float(fragment[fragment.index == "C2"].atom_y),
                  float(fragment[fragment.index == "C2"].atom_z)])

    c = np.array([float(fragment[fragment.index == "C4"].atom_x),
                  float(fragment[fragment.index == "C4"].atom_y),
                  float(fragment[fragment.index == "C4"].atom_z)])

    alpha = np.radians(109.6)

    ab, bc = b - a, c - b

    gamma = np.arccos(np.dot(ab, bc.T) / (np.linalg.norm(ab) * np.linalg.norm(bc)))

    d_angle = np.radians(180) - alpha + gamma
    cd_norm = 0.97

    # rotate with an axis perpendicular to both ab and bc
    cd = rotation_from_axis_and_angle(axis=np.cross(ab, bc), angle=d_angle, rot_vec=ab) * cd_norm

    # the first point is created, now rotate it 20 degrees each time
    dihedrals = np.arange(0, 360, 20)

    frames = []
    for i, angle in enumerate(dihedrals):
        new_point = rotation_from_axis_and_angle(axis=bc, angle=np.radians(angle), rot_vec=cd)

        # translate new point
        new_point = np.add(np.add(new_point, ab), bc)

        indexname = 'aH' + str(i + 2)

        frame = pd.DataFrame(index=[indexname],
                             data=[['H', new_point[0], new_point[1], new_point[2], settings.get_vdw_radius('H')]],
                             columns=['atom_symbol', 'atom_x', 'atom_y', 'atom_z', 'vdw_radius'])

        frames.append(copy.deepcopy(frame))

    frames.append(fragment)
    fragment = pd.concat(frames)

    return fragment


def rotation_from_axis_and_angle(axis, angle, rot_vec):
    axis = axis / np.linalg.norm(axis)
    rot_vec = rot_vec / np.linalg.norm(rot_vec)

    # the rotation is around the axis that is given, called u for now
    ux, uy, uz = axis[0], axis[1], axis[2]
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    angle = np.radians(angle)

    rot_mat = np.array([[cos_a + ux**2 * (1 - cos_a),
                         ux * uy * (1 - cos_a) - uz * sin_a,
                         ux * uz * (1 - cos_a) + uy * sin_a],
                        [uy * ux * (1 - cos_a) + uz * sin_a,
                         cos_a + uy**2 * (1 - cos_a),
                         uy * uz * (1 - cos_a) - ux * sin_a],
                        [uz * ux * (1 - cos_a) - uy * sin_a,
                         uz * uy * (1 - cos_a) + ux * sin_a,
                         cos_a + uz**2 * (1 - cos_a)]])

    return np.dot(rot_mat, rot_vec)


def average_fragment(df, settings):
    """ Returns a fragment containing the average points of the central groups. """

    settings.alignment_labels()
    central_group_df = df[df.in_central_group]

    central_group_df = central_group_df.drop(columns=['entry_id', 'id', 'in_central_group'])

    # TODO: ruthenium? if there's an R, take the average vdw
    if settings.alignment["R"] is not None:
        counts = central_group_df[central_group_df['atom_label'].str.contains("R")]['atom_symbol'].value_counts()

        vdw = 0
        atoms = 0
        for i, count in counts.items():
            atoms += count
            vdw += count * settings.get_vdw_radius(i)
        avg_vdw = vdw / atoms

    avg_fragment_df = central_group_df.groupby('atom_label').agg({'atom_symbol': 'first', 'atom_x': 'mean',
                                                                  'atom_y': 'mean', 'atom_z': 'mean'}).reset_index()

    avg_fragment_df["vdw_radius"] = 0

    for idx, row in avg_fragment_df.iterrows():
        if "R" in row.atom_label:
            avg_fragment_df.loc[idx, "vdw_radius"] = avg_vdw
        else:
            avg_fragment_df.loc[idx, "vdw_radius"] = settings.get_vdw_radius(row.atom_symbol)

    return avg_fragment_df


def get_vdw_distance_contact(df, settings):
    if settings.to_count_contact == "centroid":
        return calculate_longest_vdw_radius_contact(df, settings)

    # else return vdw radius of the atom the user is interested in
    return settings.get_vdw_radius(settings.to_count_contact)


def calculate_longest_vdw_radius_contact(df, settings):
    # TODO: if there's an R, take the biggest vdw radius
    longest_distance = 0
    atom_a = None

    # take the first fragment and it's centroid
    first_fragment_df = df[df.id == df.id.unique()[0]]
    centroid = first_fragment_df.groupby('id').mean()

    for _, atom in first_fragment_df.iterrows():
        if not atom.in_central_group:
            distance = math.sqrt((atom.atom_x - centroid.atom_x)**2 + (atom.atom_y - centroid.atom_y)**2 +
                                 (atom.atom_z - centroid.atom_z)**2)

            if distance > longest_distance:
                longest_distance = distance
                atom_a = atom

    longest_vdw_distance = (longest_distance + settings.get_vdw_radius(atom_a.atom_symbol))

    return longest_vdw_distance
