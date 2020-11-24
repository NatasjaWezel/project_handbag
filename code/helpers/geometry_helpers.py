import math

import numpy as np
import pandas as pd

import copy

from numba import jit


def make_coordinate_df(df, settings, avg_fragment):
    try:
        coordinate_df = pd.read_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

        return coordinate_df
    except (KeyError, FileNotFoundError):
        first_fragment_df = df[df.id == df.id.unique()[0]]
        find_closest_contact_atom = True

        if settings.to_count_contact == "centroid":
            # plot centroids of all contact fragments
            longest_vdw = get_vdw_distance_contact(df, settings)
            coordinate_df = df.groupby("id").mean().reset_index()

        elif len(first_fragment_df[first_fragment_df["atom_symbol"] == settings.to_count_contact]) == 1:
            longest_vdw = get_vdw_distance_contact(df, settings)

            # atom is unique, plot all of them
            coordinate_df = df[df.atom_symbol == settings.to_count_contact].reset_index().copy()

        else:
            longest_vdw = get_vdw_distance_contact(df, settings)
            coordinate_df = df[df.atom_symbol == settings.to_count_contact].reset_index().copy()

            # atom is not unique, find closest later
            find_closest_contact_atom = True

        coordinate_df = distances_closest_vdw_central(coordinate_df, avg_fragment, settings)

        if find_closest_contact_atom:
            # find closest atom
            coordinate_df = coordinate_df.loc[coordinate_df.groupby('id').distance.idxmin()].reset_index(drop=True)

        coordinate_df['longest_vdw'] = longest_vdw
        coordinate_df.to_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

        print("Done")

        return coordinate_df


def distances_closest_vdw_central(coordinate_df, avg_fragment, settings):
    length = len(coordinate_df)

    closest_distances = np.zeros(length)
    closest_atoms_vdw = np.zeros(length)

    points_avg_f = np.array([avg_fragment.atom_x, avg_fragment.atom_y, avg_fragment.atom_z]).T
    vdw_radii = np.array(avg_fragment.vdw_radius)

    print("Searching for nearest atom from contact group...")
    print(length)

    xcoord = np.array(coordinate_df.atom_x)
    ycoord = np.array(coordinate_df.atom_y)
    zcoord = np.array(coordinate_df.atom_z)

    closest_atoms_vdw, closest_distances = p_dist_calc(closest_atoms_vdw, closest_distances,
                                                       xcoord, ycoord, zcoord,
                                                       length, points_avg_f, vdw_radii)

    coordinate_df.loc[:, "distance"] = closest_distances
    coordinate_df.loc[:, "vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df


@jit(nopython=True)
def p_dist_calc(closest_atoms_vdw, closest_distances, xcoord, ycoord, zcoord, length, points_avg_f, vdw_radii):

    # use prange from numba
    for idx in range(length):

        # grab x, y and z of current contact from np arrays
        contact_point = np.array([xcoord[idx], ycoord[idx], zcoord[idx]])

        # set distance to infinite so you'll find a lower distance soon
        min_dist = 1000000000
        min_atom_vdw = None

        # calc distance with every avg fragment point, remember shortest one
        for i, avg_fragment_point in enumerate(points_avg_f):
            t_dist = np.sqrt(np.sum((avg_fragment_point - contact_point)**2, axis=0))

            # also remember the vdw radius of the closest atom
            if t_dist < min_dist:
                min_dist = t_dist
                min_atom_vdw = vdw_radii[i]

        closest_distances[idx] = min_dist
        closest_atoms_vdw[idx] = min_atom_vdw

    return closest_atoms_vdw, closest_distances


def get_dihedral_and_h(settings):
    methyl_model = {}

    with open("data/methylmodel.csv", 'r') as model_file:
        lines = model_file.readlines()

    keys = lines[0].split(",")
    labels = lines[1].split(",")

    for key, label in zip(keys, labels):
        methyl_model[key.strip()] = label.strip()

    return methyl_model


def add_model_methyl(fragment, settings):
    print("Adding model CH3 group...")

    methyl_model = get_dihedral_and_h(settings)

    h1, h2, h3 = methyl_model['h1'], methyl_model['h2'], methyl_model['h3']
    fragment = fragment[(fragment.lablabel != h1) & (fragment.lablabel != h2) & (fragment.lablabel != h3)]

    a = np.array([float(fragment[fragment.lablabel == methyl_model['dihedral1']].atom_x),
                  float(fragment[fragment.lablabel == methyl_model['dihedral1']].atom_y),
                  float(fragment[fragment.lablabel == methyl_model['dihedral1']].atom_z)])

    b = np.array([float(fragment[fragment.lablabel == methyl_model['dihedral2']].atom_x),
                  float(fragment[fragment.lablabel == methyl_model['dihedral2']].atom_y),
                  float(fragment[fragment.lablabel == methyl_model['dihedral2']].atom_z)])

    c = np.array([float(fragment[fragment.lablabel == methyl_model['dihedral3']].atom_x),
                  float(fragment[fragment.lablabel == methyl_model['dihedral3']].atom_y),
                  float(fragment[fragment.lablabel == methyl_model['dihedral3']].atom_z)])

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

        indexname = 'aH' + str(i + 1)

        frame = pd.DataFrame(data=[['H', new_point[0], new_point[1], new_point[2], settings.get_vdw_radius('H'),
                                    indexname, '-']],
                             columns=['atom_symbol', 'atom_x', 'atom_y', 'atom_z', 'vdw_radius', 'atom_label',
                                      'lablabel'])

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

    avg_fragment_df = central_group_df.groupby('atom_label').agg({'atom_symbol': 'first', 'lablabel': 'first',
                                                                  'atom_x': 'mean',
                                                                  'atom_y': 'mean',
                                                                  'atom_z': 'mean'}).reset_index()

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
    # TODO: if there's an R, kick that one out
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
