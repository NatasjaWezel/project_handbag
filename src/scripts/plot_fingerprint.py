# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `plot_fingerprint.py`
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys
import time

import pandas as pd

from classes.Settings import Settings
from classes.Radii import Radii
from classes.Fingerprint import Fingerprint

from helpers.geometry_helpers import (make_coordinate_df)
from helpers.geometry_helpers import distances_closest_vdw_central

from constants.paths import WORKDIR
from constants.constants import STANDARD_EXTRA_VDW


def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_fingerprint.py <path/to/inputfile> <reference point contact group>")
        sys.exit(1)

    t0 = time.time()

    settings = Settings(WORKDIR, sys.argv[1])
    settings.set_contact_reference_point(sys.argv[2])

    try:
        df = pd.read_csv(settings.get_aligned_csv_filename(), header=0)
        avg_frag = pd.read_csv(settings.outputfile_prefix + "_avg_fragment.csv", header=0)
    except FileNotFoundError:
        print('First align and calculate average fragment.')
        sys.exit(2)

    make_fingerprint_plots(df, avg_frag, settings, STANDARD_EXTRA_VDW)
    t1 = time.time() - t0
    print("Duration: %.2f s." % t1)


def make_fingerprint_plots(df, avg_frag, settings, STANDARD_EXTRA_VDW):
    fingerprint = Fingerprint(settings)

    radii = Radii(settings.get_radii_csv_name())
    coordinate_df = make_coordinate_df(df, settings, avg_frag, radii)
    coordinate_df['moved'] = coordinate_df['distance'] - coordinate_df['vdw_closest_atom']\
        - coordinate_df['longest_vdw']

    # make first the fingerprint plot with everything
    fingerprint.make_plot(coordinate_df, STANDARD_EXTRA_VDW)
    fingerprint.next()

    while fingerprint.not_done():
        avg_frag_f = avg_frag[avg_frag.label.isin(fingerprint.get_label_list())]
        labels = fingerprint.get_labels()

        coordinate_df_f = distances_closest_vdw_central(coordinate_df, avg_frag_f, labels)

        coordinate_df_f['moved'] = coordinate_df_f['distance' + labels]\
            - coordinate_df_f['vdw_closest_atom' + labels]\
            - coordinate_df_f['longest_vdw']

        fingerprint.make_plot(coordinate_df_f)

        fingerprint.next()


if __name__ == "__main__":
    main()
