# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script loads the coordinates of the fragments exported from a conquest query
# and aligns the central groups with the kabsch algorithm.
# It then saves the new coordinates in a .csv file.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys
import os
import time

import numpy as np
from classes.Settings import AlignmentSettings
from constants.paths import CENTRAL_GROUPS_CSV, WORKDIR
from constants.constants import MAX_DATAFILE_SIZE
from tqdm import tqdm

from helpers.general_helpers import split

from helpers.alignment_helpers import (calc_rmse, kabsch_align,
                                       perform_rotations, perform_translation,
                                       read_raw_data)


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

    split_file_if_too_big(settings.coordinate_file, settings.no_atoms)
    settings.update_coordinate_filename()

    # TODO: build check for existing alignment
    align_all_fragments(settings)

    t1 = time.time() - t0

    print("Duration: %.2f s." % t1)


def align_all_fragments(settings, to_mirror=True):
    # get filenames
    aligned_csv_filename, structures_csv_filename = settings.get_aligned_csv_filenames()

    # check if already aligned
    if os.path.exists(aligned_csv_filename):
        print("The fragments are already aligned")
        return

    # TODO: count files
    data, structures, data_matrix, first_fragment = rotate_first_fragment(settings, to_mirror)

    # align all fragments
    structures, data_matrix = do_kabsch_align(settings, data_matrix, structures, first_fragment, to_mirror)

    # put back into df
    data_matrix = data_matrix.T
    x_vec, y_vec, z_vec = data_matrix[0], data_matrix[1], data_matrix[2]
    data.x, data.y, data.z = x_vec, y_vec, z_vec

    # reindex the data to a more readable format
    data = data.reindex(['fragment_id', '_id', 'symbol', 'label', 'x', 'y', 'z'], axis=1)

    # save as csv
    structures.to_csv(structures_csv_filename, index=False)
    data.to_csv(aligned_csv_filename, index=False)


def split_file_if_too_big(filename, no_atoms):
    # do a filesize check
    filesize = os.path.getsize(filename)

    # if file bigger than max size, split it into several files
    if filesize > MAX_DATAFILE_SIZE:
        print("Splitting original datafile...")
        output_name_template = filename.rsplit('.', 1)[0] + '_%s.' + filename.rsplit('.', 1)[1]

        # split after fragments, each fragment has no_atoms + 1 header line
        row_limit = 6e6 - (6e6 % (no_atoms + 1))

        # TODO: make sure it doesn't split if file is small enough
        split(open(filename), delimiter=',', row_limit=row_limit,
              output_name_template=output_name_template, output_path='.')


def rotate_first_fragment(settings, to_mirror):
    data, structures = read_raw_data(settings.coordinate_file, settings.no_atoms)

    # bin atoms to be binned
    settings, fragments, data = prepare_data(settings, data)

    # restructure df to matrix
    data_matrix = np.array([np.array(data.x), np.array(data.y), np.array(data.z)]).T
    A = data_matrix[:settings.no_atoms]

    # translate and rotate first fragment onto the origin as for nice viewing
    A = perform_translation(A.copy(), settings.get_index_alignment_atom('center'))
    A = perform_rotations(A.copy(), [settings.get_index_alignment_atom('yaxis'),
                                     settings.get_index_alignment_atom('xyplane')])

    if to_mirror:
        mirrored, A = mirror(A)

        if mirrored:
            structures.loc[structures.index == 0, 'mirrored'] = True

    data_matrix[:settings.no_atoms] = A

    A = data_matrix[0:settings.no_atoms_central]

    return data, structures, data_matrix, A


def mirror(matrix):
    # if mean of z is negative, reflect the whole fragment
    if matrix[:, 2].mean() < 0:
        # mirror by switching signs of z coordinate
        matrix[:, 2] *= -1

        return True, matrix

    return False, matrix


def do_kabsch_align(settings, data_matrix, structures, A, to_mirror):
    no_atoms = settings.no_atoms
    no_atoms_central = settings.no_atoms_central

    n = A.shape[0]

    print("Applying Kabsch Algorithm...")
    for i in tqdm(range(1, settings.get_no_fragments())):
        B_central = data_matrix[i * no_atoms:no_atoms_central + i * no_atoms]
        B_total = data_matrix[i * no_atoms: (i + 1) * no_atoms]

        # run kabsch, shift frame in data each time
        B_total_2 = kabsch_align(A, B_central, B_total, n)

        if to_mirror:
            mirrored, B_total_2 = mirror(B_total_2)

            if mirrored:
                structures.loc[structures.index == i, 'mirrored'] = True

        # calculate and save error
        rmse = calc_rmse(A, B_total_2[:no_atoms_central], n)
        structures.loc[structures.index == i, 'rmse'] = rmse

        # put back into overall matrix
        data_matrix[i * no_atoms: (i + 1) * no_atoms] = B_total_2

    return structures, data_matrix


def prepare_data(settings, data):
    amount_rows = len(data)

    # calc amount of fragments
    no_fragments = int(amount_rows / settings.no_atoms)

    settings.set_no_fragments(no_fragments)

    # give each fragment a unique id
    fragment_ids = range(0, no_fragments)
    fragment_ids = np.repeat(fragment_ids, settings.no_atoms)
    data['fragment_id'] = fragment_ids

    # give each atom in each fragment their label from conquest
    labels = settings.label_list * no_fragments
    data['label'] = labels

    if settings.alignment['bin'] != '-':
        print(f"Throwing away {len(settings.alignment['bin'])} atoms.")
        settings.no_atoms -= len(settings.alignment['bin'])
        settings.no_atoms_central -= len(settings.alignment['bin'])

        for label in settings.alignment['bin']:
            settings.label_list.remove(label)

        data = data[~data.label.isin(settings.alignment['bin'])].reset_index()

    return settings, no_fragments, data


if __name__ == "__main__":
    main()
