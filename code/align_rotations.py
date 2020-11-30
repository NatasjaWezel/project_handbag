# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# It loads the coordinates of the fragments exported from a conquest query and
# aligns the central groups by using rotation matrices and other linear algebra.
# It then saves the new coordinates in a .csv file.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys
import time

import numpy as np
from tqdm import tqdm

from helpers.alignment_helpers import (alignment_dict, calc_rmse,
                                       perform_rotations, perform_translation,
                                       read_coord_file, read_raw_data)
from align_kabsch import prepare_data


def main():

    if len(sys.argv) != 2:
        print("Usage: python load_from_coords.py <path/to/inputfile>")
        sys.exit(1)

    t0 = time.time()

    filename = sys.argv[1]

    title = filename.rsplit('\\')[-1].rsplit('.', 1)[0]
    central = title.split("_")[0]
    contact = title.split("_")[1]

    # read part of the datafile first for preparation
    no_atoms, no_atoms_central, label_list = read_coord_file(filename)

    data, structures = read_raw_data(filename, no_atoms)
    fragments, data = prepare_data(data, no_atoms, label_list)
    
    alignment = alignment_dict(central_group_name=filename.rsplit('\\')[-1].rsplit('.', 1)[0]
                               .rsplit('_aligned', 1)[0].split("_")[0])

    # restructure df to matrix
    data_xyz_matrix = np.array([np.array(data.x), np.array(data.y), np.array(data.z)]).T

    # rotate first fragment
    first_matrix = data_xyz_matrix[0:no_atoms]
    first_matrix = perform_translation(first_matrix, index_center=label_list.index(alignment['center']))
    first_matrix = perform_rotations(first_matrix.copy(),
                                                    [label_list.index(alignment['yaxis']),
                                                    label_list.index(alignment['xyplane'])])

    if first_matrix[:,2].mean() < 0:
        # mirror by switching signs of z coordinate
        first_matrix[:,2] *= -1
        structures.loc[structures.index == 0, 'mirrored'] = True

    # put back into the matrix
    data_xyz_matrix[0:no_atoms] = first_matrix

    # to calculate rmse we only need the central group
    first_matrix = first_matrix[:no_atoms]
    n = first_matrix.shape[0]

    # ***second part***: all the other fragments
    print("Rotating all fragments...")
    for i in tqdm(range(1, fragments)):
        # translate and rotate first fragment onto the origin as for nice viewing
        begin, end = i * no_atoms, (i + 1) * no_atoms
        part_matrix = data_xyz_matrix[begin:end].copy()

        # give the index of which atom is needed for the calculations
        part_matrix = perform_translation(part_matrix, index_center=label_list.index(alignment['center']))

        part_matrix = perform_rotations(part_matrix.copy(),
                                                       [label_list.index(alignment['yaxis']),
                                                        label_list.index(alignment['xyplane'])])

        # invert if necessary
        if part_matrix[:,2].mean() < 0:
            # mirror by switching signs of z coordinate
            part_matrix[:,2] *= -1
            structures.loc[structures.index == i, 'mirrored'] = True

        data_xyz_matrix[begin:end] = part_matrix

        # calculate and save error
        rmse = calc_rmse(first_matrix, part_matrix[:no_atoms], n)
        structures.loc[structures.index == i, 'rmse',] = rmse

    # put back into df
    data_xyz_matrix_T = data_xyz_matrix.T
    x_vec, y_vec, z_vec = data_xyz_matrix_T[0], data_xyz_matrix_T[1], data_xyz_matrix_T[2]
    data.x, data.y, data.z = x_vec, y_vec, z_vec

    # reindex the data to a more readable format
    data = data.reindex(['fragment_id', '_id', 'symbol', 'label', 'x', 'y', 'z'], axis=1)

    # save rmse's to csv
    structures.to_csv('results/' + central + '/' + central + "_" + contact + "_vdw.5/"\
                      + central + "_" + contact + "_rotation_structures_12-4-3.csv", index=False)

    # save aligned coordinates to csv
    data.to_csv('results/' + central + '/' + central + "_" + contact + "_vdw.5/"\
                + central + "_" + contact + "_rotation_aligned_12-4-3.csv", index=False)

    t1 = time.time() - t0

    print("Duration: %.2f s." % t1)


if __name__ == "__main__":
    main()
