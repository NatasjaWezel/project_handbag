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


import pandas as pd

from classes.Settings import Settings
from classes.Radii import Radii
from helpers.density_helpers import make_density_df, find_available_volume
from helpers.geometry_helpers import make_coordinate_df

from constants.paths import WORKDIR


def main():

    if len(sys.argv) != 4:
        print("Usage: python plot_density.py <path/to/inputfile> <resolution> <contact group reference point>")
        sys.exit(1)

    t0 = time.time()

    settings = Settings(WORKDIR, sys.argv[1])
    settings.set_contact_reference_point(sys.argv[3])

    # resolution of the bins, in Angstrom
    settings.set_resolution(float(sys.argv[2]))

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
    tolerance = 0.5
    contact_group_radius = radii.get_vdw_distance_contact(settings.contact_rp)
    Vavailable = find_available_volume(avg_fragment=avg_frag, extra=(tolerance + contact_group_radius))
    print('Available volume:', Vavailable)

    threshold_calc = density_df.datafrac_normalized.max() * settings.threshold

    in_cluster = density_df[density_df.datafrac_normalized >= threshold_calc]
    datafrac = in_cluster.datafrac_normalized.sum()

    Vcluster = len(in_cluster) * settings.resolution**3

    directionality = datafrac / Vcluster * (Vavailable/2)

    print(f"Directionality: {directionality}")
    print(density_df[settings.contact_rp].sum(), datafrac, Vcluster, Vavailable)


if __name__ == "__main__":
    main()
