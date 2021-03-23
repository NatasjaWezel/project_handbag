# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `Radii` is a class that takes as input the Radii csv, and is used to get the covalent and vanderwaals radii of
# elements. They are saved in a dictionary, so that if the same element is looked up again, we do not need to open the
# file and read it again.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pandas as pd


class Radii():
    def __init__(self, RADII_CSV):
        self.radii_filename = RADII_CSV
        self.vdw_radii = {}
        self.cov_radii = {}

    def get_vdw_radius(self, symbol):
        """ Adds element with vdw radius to dictionary if not exists yet, and returns its vdw radius. """

        # if not looked up before
        if symbol not in self.vdw_radii.keys():

            # read file, find and add vdw_radius to dictionary
            radii_df = pd.read_csv(self.radii_filename, comment="#")
            test_df = radii_df[radii_df.symbol == symbol]

            assert len(test_df) == 1, "You're trying to look up vdw radius of element that is not in bondi's list."

            self.vdw_radii[symbol] = float(test_df.vdw_radius)

        # return dictionary value
        return self.vdw_radii[symbol]

    def get_cov_radius(self, symbol):
        """ Adds element with cov radius to dictionary if not exist yet, and returns its covalent radius. """

        # if not looked up before
        if symbol not in self.cov_radii.keys():

            # read file, find and add cov_radius to dictionary
            radii_df = pd.read_csv(self.radii_filename, comment="#")
            test_df = radii_df[radii_df.symbol == symbol]

            assert len(test_df) == 1, "You're trying to look up cov radius of element that is not in bondi's list."

            self.cov_radii[symbol] = float(test_df.cov_radius)

        # return dictionary value
        return self.cov_radii[symbol]

    def get_vdw_distance_contact(self, contact_rp):
        """ Returns the vanderwaals radius from the atom that is the reference point of the contact group. If the rp is
            the centroid of a ring, returns the carbon vdw radius. """

        if contact_rp.lower() == "centroid":
            return self.get_vdw_radius("C")

        # else return vdw radius of the atom the user is interested in
        return self.get_vdw_radius(contact_rp)
