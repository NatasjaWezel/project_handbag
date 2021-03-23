# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `Fingerprint` is a class that takes as input the central group, reads the fingerprint.csv and can generate all the
# fingerprint plots by using its 'next' function.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
import os


class Fingerprint():
    """ This class is used to make all fingerprint plots for a given central group, provided by the input argument
        settings. It uses the information in the fingerprints.csv to see what labels to use. """

    def __init__(self, settings):
        """ Initialize the fingerprint plot with settings and check if all files exist. """

        self.central = settings.central_name
        self.contact = settings.contact_name
        self.contact_rp = settings.contact_rp

        # get csv filename from settings object
        self.csv = settings.get_fingerprint_filename()
        self.WORKDIR = settings.WORKDIR

        # initialize with -1 so the first plot can be made if counter == -1
        self.counter = -1

        # read the csv and get only specific information for this central group
        df = pd.read_csv(self.csv, comment="#")
        self.specific = df[df.central == self.central].reset_index()

        # check if the result file folders
        if not os.path.exists(self.WORKDIR + '/results/'):
            os.mkdir(self.WORKDIR + '/results/')
            os.mkdir(self.WORKDIR + '/results/pairs')
            os.mkdir(f'{self.WORKDIR}/results/pairs/{self.central}/')
        elif not os.path.exists(self.WORKDIR + '/results/pairs'):
            os.mkdir(self.WORKDIR + '/results/pairs')
            os.mkdir(f'{self.WORKDIR}/results/pairs/{self.central}/')
        elif not os.path.exists(f'{self.WORKDIR}/results/pairs/{self.central}/'):
            os.mkdir(f'{self.WORKDIR}/results/pairs/{self.central}/')

    def get_labels(self):
        """ Get labels from file, except if the first plot is made. """

        if self.counter == -1:
            return ""

        return self.specific[self.specific.index == self.counter].labels.item()

    def get_label_list(self):
        """ If list contains multiple labels, get list."""

        return self.specific[self.specific.index == self.counter].labels.item().split('&')

    def get_description(self):
        """ Get description from file, except is the first plot is made. """

        if self.counter == -1:
            return "closest atom"

        return self.specific.iloc[self.counter].description

    def not_done(self):
        """ A boolean to see if all plots are made or not. """

        return self.counter < len(self.specific)

    def next(self):
        """ Increment plot counter if next plot is being made. """

        self.counter += 1

    def make_plot(self, coordinate_df, STANDARD_EXTRA_VDW):
        """ Makes the plots. Takes as input the coordinate df with the vdw corrected distance already calculated. """

        fig = plt.figure(figsize=(8, 4))
        fig.subplots_adjust(bottom=0.3)
        plt.title(f"Fingerprint of {self.central} ({self.get_description()})--{self.contact} ({self.contact_rp})")

        # get two dataframes: one containing the points with vdw overlap and one without
        test_neg = coordinate_df[coordinate_df["moved"] < 0]
        test_pos = coordinate_df[coordinate_df["moved"] >= 0]

        # write down information in the plot
        plt.figtext(0.15, 0.11, f"Negative fraction: {len(test_neg)/len(coordinate_df) * 100 :.2f}%,\
                    Mean: {test_neg['moved'].mean() :.2f}$\\AA$")
        plt.figtext(0.15, 0.06, f"Positive fraction: {len(test_pos)/len(coordinate_df) * 100 :.2f}%,\
                    Mean: {test_pos['moved'].mean() :.2f}$\\AA$")
        plt.figtext(0.15, 0.01, f"Overall mean: {coordinate_df['moved'].mean() :.2f}$\\AA$")

        plt.xlabel("VDW overlap ($\\AA$)")
        plt.ylabel("Percentage of total data points")

        plt.grid('on')

        # set xlim for easy comparing, change this if this range is not suitable
        plt.xlim(-2, 3)

        # give the lines for vdw corrected distance, 0 and the maximum distance
        plt.vlines(0, 0, 0.15, color="black", label="VDW radius atom central")
        plt.vlines(coordinate_df['vdw_closest_atom'].max(), 0, 0.15, color="lightgreen", label="VDW radii")
        plt.vlines(coordinate_df['vdw_closest_atom'].max() + STANDARD_EXTRA_VDW, 0, 0.15,
                   color="green",
                   label="VDW radii + " + str(STANDARD_EXTRA_VDW))

        # calculate and normalize the heights of the bins
        heights, bins = np.histogram(coordinate_df.moved, bins='auto')
        heights = (heights/sum(heights) * 100)

        plt.bar(bins[:-1], heights, width=(max(bins) - min(bins))/len(bins)+0.01)

        plt.legend()

        title = f'{self.WORKDIR}/results/pairs/{self.central}/{self.central}_{self.contact}_{self.contact_rp}'
        title += f'_fingerprint_{self.get_labels()}.png'
        plt.savefig(title)
        plt.close(fig)

        print(f"Fingerprint saved in {title}")
