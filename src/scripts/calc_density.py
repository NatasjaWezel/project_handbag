# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments. It then divides the
# surrounding space into a number of bins, depending on which resolution is
# set. It counts how many of the contact atoms/ centers of contact groups are
# are in each bin and normalizes that by the total amount of contact atoms or
# groups.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys
import time

import numpy as np
import pandas as pd

from classes.Settings import Settings
from classes.Radii import Radii
from helpers.density_helpers import make_density_df, find_available_volume, calc_distances
from helpers.geometry_helpers import make_coordinate_df

from constants.constants import STANDARD_THRESHOLD, STANDARD_RES, STANDARD_EXTRA_VDW
from constants.paths import WORKDIR


def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_density.py <path/to/inputfile> <contact group reference point>")
        sys.exit(1)

    t0 = time.time()

    settings = Settings(WORKDIR, sys.argv[1])
    settings.set_contact_reference_point(sys.argv[2])

    # resolution of the bins, in Angstrom
    settings.set_resolution(STANDARD_RES)
    settings.set_threshold(STANDARD_THRESHOLD)

    try:
        df = pd.read_csv(settings.get_aligned_csv_filename(), header=0)
        avg_frag = pd.read_csv(settings.outputfile_prefix + "_avg_fragment.csv", header=0)
    except FileNotFoundError:
        print('First align and calculate average fragment.')
        sys.exit(2)

    radii = Radii(settings.get_radii_csv_name())
    # grab only the atoms that are in the contact groups
    coordinate_df = make_coordinate_df(df, settings, avg_frag, radii)

    density_df = make_density_df(settings, coordinate_df, again=True)

    t1 = time.time() - t0
    print("Duration: %.2f s." % t1)

    # find the volume of the central group
    tolerance = STANDARD_EXTRA_VDW
    contact_group_radius = radii.get_vdw_distance_contact(settings.contact_rp)
    Vavailable = find_available_volume(avg_fragment=avg_frag, extra=(tolerance + contact_group_radius))
    print('Available volume:', Vavailable)

    threshold_calc = density_df.datafrac_normalized.max() * settings.threshold

    in_cluster = density_df[density_df.datafrac_normalized >= threshold_calc]
    datafrac = in_cluster.datafrac_normalized.sum()

    Vcluster = len(in_cluster) * settings.resolution**3

    directionality = datafrac / Vcluster * (Vavailable/2)

    print(f"Directionality: {directionality}")

    vdw_overlap = calc_vdw_overlap(in_cluster.copy(), settings, avg_frag, contact_group_radius)
    print(f"Vdw overlap datapoints in cluster: {vdw_overlap :.2f}%")


def calc_vdw_overlap(in_cluster, settings, avg_fragment, contact_group_radius):

    # center of bin is at xstart + a half times the resolution
    in_cluster['x_center'] = in_cluster.xstart + 0.5 * settings.resolution
    in_cluster['y_center'] = in_cluster.ystart + 0.5 * settings.resolution
    in_cluster['z_center'] = in_cluster.zstart + 0.5 * settings.resolution

    bin_coordinates = np.transpose(np.array([in_cluster.x_center, in_cluster.y_center, in_cluster.z_center]))
    in_vdw_vol = np.zeros(len(in_cluster))

    for i, atom in avg_fragment.iterrows():
        indices = np.transpose(np.where(in_vdw_vol == 0))

        fragment_point = np.array([atom.x, atom.y, atom.z, atom.vdw_radius])

        in_vdw_vol = calc_distances(in_vdw_vol, bin_coordinates, fragment_point, indices,
                                    extra=(contact_group_radius))

    in_cluster['in_vdw_vol'] = in_vdw_vol

    return in_cluster[in_cluster.in_vdw_vol != 0].datafrac_normalized.sum() * 100


if __name__ == "__main__":
    main()
