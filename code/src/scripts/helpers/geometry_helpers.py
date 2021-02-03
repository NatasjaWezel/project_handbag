import numpy as np
import pandas as pd

import copy

import time

from numba import jit


def make_coordinate_df(df, settings, avg_fragment, radii):
    try:
        coordinate_df = pd.read_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())
        print("Coordinate df already existed, loaded from file")
    except (KeyError, FileNotFoundError):
        print("Searching for nearest atom from central group...")
        t0 = time.time()

        df = df[df.label == "-"]

        first_fragment_df = df[df.fragment_id == df.fragment_id.unique()[0]]

        print("Atoms in contact group:", len(first_fragment_df), "atom to count: ", settings.to_count_contact)
        find_closest_contact_atom = False

        if settings.to_count_contact == "centroid":
            # plot centroids of all contact fragments
            longest_vdw = radii.get_vdw_radius("C")
            coordinate_df = df.groupby("fragment_id").mean().reset_index()

        elif len(first_fragment_df[first_fragment_df["symbol"] == settings.to_count_contact]) == 1:
            longest_vdw = radii.get_vdw_distance_contact(df, settings)

            # atom is unique, plot all of them
            coordinate_df = df[df.symbol == settings.to_count_contact].reset_index().copy()

        else:
            longest_vdw = radii.get_vdw_distance_contact(df, settings)
            coordinate_df = df[df.symbol == settings.to_count_contact].reset_index().copy()

            # atom is not unique, find closest later
            find_closest_contact_atom = True

        coordinate_df = distances_closest_vdw_central(coordinate_df, avg_fragment)

        if find_closest_contact_atom:
            # find closest atom
            coordinate_df = coordinate_df.loc[coordinate_df.groupby('fragment_id').distance.idxmin()]\
                                         .reset_index(drop=True)

        coordinate_df['longest_vdw'] = longest_vdw
        coordinate_df.to_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

        t1 = time.time()
        print("Coordinate df is made, duration:", t1-t0, 's')

    return coordinate_df


def distances_closest_vdw_central(coordinate_df, avg_fragment):
    length = len(coordinate_df)

    closest_distances = np.zeros(length)
    closest_atoms_vdw = np.zeros(length)

    points_avg_f = np.array([avg_fragment.x, avg_fragment.y, avg_fragment.z]).T
    vdw_radii = np.array(avg_fragment.vdw_radius)

    xcoord = np.array(coordinate_df.x)
    ycoord = np.array(coordinate_df.y)
    zcoord = np.array(coordinate_df.z)

    closest_atoms_vdw, closest_distances = p_dist_calc(closest_atoms_vdw, closest_distances,
                                                       xcoord, ycoord, zcoord,
                                                       length, points_avg_f, vdw_radii)

    coordinate_df.loc[:, "distance"] = closest_distances
    coordinate_df.loc[:, "vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df


@jit(nopython=True)
def p_dist_calc(closest_atoms_vdw, closest_distances, xcoord, ycoord, zcoord, length, points_avg_f, vdw_radii):

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


def get_dihedral_and_h(CSV, central_group_name):
    methyl_model = {}

    df = pd.read_csv(CSV, header=0)

    df = df[df.central_group == central_group_name]

    methyl_model["dihedral1"] = df.dihedral1.item()
    methyl_model["dihedral2"] = df.dihedral2.item()
    methyl_model["dihedral3"] = df.dihedral3.item()

    return methyl_model


def add_model_methyl(CSV, fragment, settings, radii):
    print("Adding model CH3 group...", end=" ")

    methyl_model = get_dihedral_and_h(CSV, settings.central_group_name)

    # if labels have been switched by kmeans, look at old labels
    column = "label"

    a = np.array([float(fragment[fragment[column].str.contains(methyl_model['dihedral1'])].x),
                  float(fragment[fragment[column].str.contains(methyl_model['dihedral1'])].y),
                  float(fragment[fragment[column].str.contains(methyl_model['dihedral1'])].z)])

    b = np.array([float(fragment[fragment[column].str.contains(methyl_model['dihedral2'])].x),
                  float(fragment[fragment[column].str.contains(methyl_model['dihedral2'])].y),
                  float(fragment[fragment[column].str.contains(methyl_model['dihedral2'])].z)])

    c = np.array([float(fragment[fragment[column].str.contains(methyl_model['dihedral3'])].x),
                  float(fragment[fragment[column].str.contains(methyl_model['dihedral3'])].y),
                  float(fragment[fragment[column].str.contains(methyl_model['dihedral3'])].z)])

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
        new_point = np.add(np.add(np.add(new_point, a), ab), bc)

        indexname = 'aH' + str(i + 1)

        frame = pd.DataFrame(data=[['H', indexname, new_point[0], new_point[1], new_point[2],
                                    radii.get_vdw_radius('H')]],
                             columns=['symbol', 'label', 'x', 'y', 'z', 'vdw_radius'])

        frames.append(copy.deepcopy(frame))

    frames.append(fragment)
    fragment = pd.concat(frames)

    print("Done")

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


def average_fragment(df, settings, radii):
    """ Returns a fragment containing the average points of the central groups. """

    central_group_df = df[df['label'] != '-']

    # take average R vdw radius
    counts = central_group_df[central_group_df['label'].str.contains("R")]['symbol'].value_counts()

    if len(counts) > 0:
        # TODO: what happens if multiple R?
        if 'kmeans_label' not in central_group_df.columns:
            print("\nR consists of:")
            elements = counts.index.to_list()
            counts_list = counts.to_list()
            percentages = [count/np.sum(counts) for count in counts_list]

            for element in elements:
                print(element.ljust(10), end="")
            print()
            for percentage in percentages:
                print(f"{percentage * 100 :.2f}%    ".ljust(10), end="")
            print('\n')

        vdw, cov = 0, 0
        atoms = 0
        for i, count in counts.items():
            atoms += count
            vdw += count * radii.get_vdw_radius(i)
            cov += count * radii.get_cov_radius(i)
        avg_vdw = vdw / atoms
        avg_cov = cov / atoms

    if 'kmeans_label' not in central_group_df.columns:
        # sort must be false to preserve the order of the rows/labels
        avg_fragment_df = central_group_df.groupby('label', sort=False).agg({'symbol': 'first',
                                                                             'x': 'mean',
                                                                             'y': 'mean',
                                                                             'z': 'mean'}).reset_index()
    else:
        # sort must be false to preserve the order of the rows/labels
        avg_fragment_df = central_group_df.groupby('kmeans_label', sort=False).agg({'symbol': 'first',
                                                                                    'label': 'first',
                                                                                    'x': 'mean',
                                                                                    'y': 'mean',
                                                                                    'z': 'mean'}).reset_index()

    avg_fragment_df["vdw_radius"] = 0
    avg_fragment_df["cov_radius"] = 0

    for idx, row in avg_fragment_df.iterrows():
        if "R" in row.label:
            avg_fragment_df.loc[idx, "vdw_radius"] = avg_vdw
            avg_fragment_df.loc[idx, "cov_radius"] = avg_cov
        else:
            avg_fragment_df.loc[idx, "vdw_radius"] = radii.get_vdw_radius(row.symbol)
            avg_fragment_df.loc[idx, "cov_radius"] = radii.get_cov_radius(row.symbol)

    return avg_fragment_df
