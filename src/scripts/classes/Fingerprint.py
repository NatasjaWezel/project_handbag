import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
import os


class Fingerprint():
    def __init__(self, settings):
        self.central = settings.central_name
        self.contact = settings.contact_name
        self.contact_rp = settings.contact_rp

        self.csv = settings.get_finger_print_filename()

        self.counter = -1

        df = pd.read_csv(self.csv)
        self.specific = df[df.central == self.central].reset_index()

        if not os.path.exists('../../results/'):
            os.mkdir('../../results')
            os.mkdir('../../results/fingerprints')
            os.mkdir(f'../../results/fingerprints/{self.central}/')
        elif not os.path.exists('../../results/fingerprints'):
            os.mkdir('../../results/fingerprints')
            os.mkdir(f'../../results/fingerprints/{self.central}/')
        elif not os.path.exists(f'../../results/fingerprints/{self.central}/'):
            os.mkdir(f'../../results/fingerprints/{self.central}/')

    def get_labels(self):
        if self.counter == -1:
            return ""

        return self.specific[self.specific.index == self.counter].labels.item()

    def get_label_list(self):
        return self.specific[self.specific.index == self.counter].labels.item().split('&')

    def get_description(self):
        if self.counter == -1:
            return "closest atom"
        return self.specific.iloc[self.counter].description

    def not_done(self):
        return self.counter < len(self.specific)

    def next(self):
        self.counter += 1

    def make_plot(self, coordinate_df):
        fig = plt.figure(figsize=(8, 4))
        fig.subplots_adjust(bottom=0.3)
        plt.title(f"Fingerprint of {self.central} ({self.get_description()})--{self.contact} ({self.contact_rp})")

        test_neg = coordinate_df[coordinate_df["moved"] < 0]
        test_pos = coordinate_df[coordinate_df["moved"] >= 0]

        plt.figtext(0.15, 0.11, f"Negative fraction: {len(test_neg)/len(coordinate_df) * 100 :.2f}%,\
                    Mean: {test_neg['moved'].mean() :.2f}$\\AA$")
        plt.figtext(0.15, 0.06, f"Positive fraction: {len(test_pos)/len(coordinate_df) * 100 :.2f}%,\
                    Mean: {test_pos['moved'].mean() :.2f}$\\AA$")
        plt.figtext(0.15, 0.01, f"Overall mean: {coordinate_df['moved'].mean() :.2f}$\\AA$")

        plt.xlabel("VDW overlap ($\\AA$)")
        plt.ylabel("Fraction")

        plt.grid('on')

        plt.xlim(-2, 3)

        plt.vlines(0, 0, 0.15, color="black", label="VDW radius atom central")

        plt.vlines(coordinate_df['vdw_closest_atom'].max(), 0, 0.15, color="lightgreen", label="VDW radii")
        plt.vlines(coordinate_df['vdw_closest_atom'].max() + 0.5, 0, 0.15, color="green", label="VDW radii + 0.5")

        heights, bins = np.histogram(coordinate_df.moved, bins='auto')
        heights = (heights/sum(heights))

        plt.bar(bins[:-1], heights, width=(max(bins) - min(bins))/len(bins)+0.01)

        plt.legend()

        title = f'../../results/fingerprints/{self.central}/{self.central}_{self.contact}_{self.contact_rp}'
        title += f'_fingerprint_{self.get_labels()}.png'
        plt.savefig(title)
        plt.close(fig)

        print(f"Fingerprint saved in {title}")
