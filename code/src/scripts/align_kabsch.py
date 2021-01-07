# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script loads the coordinates of the fragments exported from a conquest query
# and aligns the central groups with the kabsch algorithm.
# It then saves the new coordinates in a .csv file.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys

import time

import numpy as np
from classes.Settings import AlignmentSettings
from constants.paths import CENTRAL_GROUPS_CSV, WORKDIR
from tqdm import tqdm

from helpers.alignment_helpers import (alignment_dict, calc_rmse, kabsch_align,
                                       perform_rotations, perform_translation,
                                       read_coord_file, read_raw_data)


def main():

    if len(sys.argv) != 2:
        print("Usage: python load_from_coords.py <path/to/coordinatefile>")
        sys.exit(1)

    t0 = time.time()

    coordinate_file = sys.argv[1]
    labelfile = coordinate_file.rsplit('.', 1)[0] + '.csv'

    settings = AlignmentSettings(WORKDIR, coordinate_file, labelfile)

    settings.set_central_group_csv(CENTRAL_GROUPS_CSV)
    settings.prepare_alignment()

    do_kabsch_align(settings)

    t1 = time.time() - t0

    print("Duration: %.2f s." % t1)


def do_kabsch_align(settings):
    no_atoms = settings.no_atoms

    data, structures = read_raw_data(settings.coordinate_data, settings.no_atoms)
    fragments, data = prepare_data(data, settings)

    # restructure df to matrix
    data_xyz_matrix = np.array([np.array(data.x), np.array(data.y), np.array(data.z)]).T

    A = data_xyz_matrix[:no_atoms]
    # translate and rotate first fragment onto the origin as for nice viewing
    A = perform_translation(A.copy(), settings.get_index_alignment_atom('center'))
    A = perform_rotations(A.copy(), [settings.get_index_alignment_atom('center')),
                                     settings.get_index_alignment_atom('center')])

    # if mean of z is negative, reflect the whole fragment
    if A[:, 2].mean() < 0:
        # mirror by switching signs of z coordinate
        A[:, 2] *= -1
        structures.loc[structures.index == 0, 'mirrored'] = True

    data_xyz_matrix[:no_atoms] = A

    A = data_xyz_matrix[0:no_atoms_central]
    n = A.shape[0]

    print("Applying Kabsch Algorithm...")
    for i in tqdm(range(1, fragments)):
        B_central = data_xyz_matrix[i * no_atoms:no_atoms_central + i * no_atoms]
        B_total = data_xyz_matrix[i * no_atoms: (i + 1) * no_atoms]

        # run kabsch, shift frame in data each time
        B_total_2 = kabsch_align(A, B_central, B_total, n)

        # invert if necessary
        if B_total_2[:, 2].mean() < 0:
            # mirror by switching signs of z coordinate
            B_total_2[:, 2] *= -1
            structures.loc[structures.index == i, 'mirrored'] = True

        # calculate and save error
        rmse = calc_rmse(A, B_total_2[:no_atoms_central], n)
        structures.loc[structures.index == i, 'rmse'] = rmse

        # put back into overall matrix
        data_xyz_matrix[i * no_atoms: (i + 1) * no_atoms] = B_total_2

    # put back into df
    data_xyz_matrix_T = data_xyz_matrix.T
    x_vec, y_vec, z_vec = data_xyz_matrix_T[0], data_xyz_matrix_T[1], data_xyz_matrix_T[2]
    data.x, data.y, data.z = x_vec, y_vec, z_vec

    # reindex the data to a more readable format
    data = data.reindex(['fragment_id', '_id', 'symbol', 'label', 'x', 'y', 'z'], axis=1)

    aligned_csv_filename, structures_csv_filename = settings.get_aligned_csv_filenames()

    # save as csv
    structures.to_csv(structures_csv_filename, index=False)
    data.to_csv(aligned_csv_filename, index=False)


def prepare_data(settings, data):
    amount_rows = len(data)

    # calc amount of fragments
    fragments = int(amount_rows / settings.no_atoms)

    # give each fragment a unique id
    fragment_ids = range(0, fragments)
    fragment_ids = np.repeat(fragment_ids, settings.no_atoms)
    data['fragment_id'] = fragment_ids

    # give each atom in each fragment their label from conquest
    labels = settings.label_list * fragments
    data['label'] = labels

    return fragments, data


if __name__ == "__main__":
    main()
