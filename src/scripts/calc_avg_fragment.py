import sys
import pandas as pd

from classes.Settings import Settings
from classes.Radii import Radii
from helpers.geometry_helpers import average_fragment, add_model_methyl
from helpers.alignment_helpers import calc_rmse

from constants.constants import RMSE_TEST

from sklearn.cluster import KMeans

from constants.paths import WORKDIR

import numpy as np


def main():

    if len(sys.argv) != 2:
        print("Usage: python calc_avg_fragment_2.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    avg_frag_settings = Settings(WORKDIR, inputfilename)
    df = pd.read_csv(avg_frag_settings.get_aligned_csv_filename(), header=0)

    # make radii object to get vdw radii
    radii = Radii(avg_frag_settings.get_radii_csv_name())

    fragment = calc_avg_frag(df, avg_frag_settings, radii)

    # get name and save average fragment
    avg_frag_file = avg_frag_settings.get_avg_frag_filename()
    fragment.to_csv(avg_frag_file, index=False)


def calc_avg_frag(df, avg_frag_settings, radii):
    fragment = average_fragment(df, avg_frag_settings, radii)

    # test
    calc_kabsch_rmse(avg_frag_settings)
    rmse_avg = calc_avg_rmse(fragment, avg_frag_settings)

    # dependent on rmses, do kmeans or not
    if rmse_avg > RMSE_TEST:
        print("RMSEs too high. Resetting labels using KMeans")
        df = reset_labels_with_kmeans(df, avg_frag_settings)

        # sort df based on new labels
        df = df.sort_values(['fragment_id', 'kmeans_label'])

        fragment = average_fragment(df, avg_frag_settings, radii)

        rmse_avg = calc_avg_rmse(fragment, avg_frag_settings, df)

    if avg_frag_settings.central_name in list(pd.read_csv(avg_frag_settings.get_methyl_csv_filename(),
                                              header=0).central):
        fragment = add_model_methyl(CSV=avg_frag_settings.get_methyl_csv_filename(), fragment=fragment,
                                    settings=avg_frag_settings, radii=radii)

    return fragment


def reset_labels_with_kmeans(df, avg_frag_settings):
    # only keep central group
    df = df[df.label != "-"].copy().reset_index()

    first_frag = df[df.fragment_id == 0]
    no_atoms_central = len(first_frag)

    first_frag = np.array([np.array(first_frag.x), np.array(first_frag.y), np.array(first_frag.z)]).T
    data = np.array([np.array(df.x), np.array(df.y), np.array(df.z)]).T

    # initialize kmeans
    kmeans = KMeans(n_clusters=no_atoms_central, init=first_frag, n_init=1)

    kmeans.fit(data)
    labels = kmeans.labels_

    df['kmeans_label'] = [str(label) for label in labels]

    return df


def calc_kabsch_rmse(settings):
    structures_file = settings.get_structure_csv_filename()

    first_rmse = pd.read_csv(structures_file).rmse.mean()
    print(f"Average RMSE kabsch alignment: {first_rmse :.2f}")

    return first_rmse


def get_aligned_fragments_for_rmse(settings):
    frags_file = settings.get_aligned_csv_filename()

    # this is the old number of atoms
    rows = 5000

    df = pd.read_csv(frags_file, nrows=rows)

    # correct for atoms that are thrown away after reading
    no_atoms = len(df[df.fragment_id == 0])
    rows = 100 * no_atoms
    df = pd.read_csv(frags_file, nrows=rows)

    return df


def calc_avg_rmse(avg_fragment, settings, df=None):

    if df is None:
        df = get_aligned_fragments_for_rmse(settings)

    no_atoms_central = len(df[(df.fragment_id == 0) & (df.label != "-")])

    unique_frags = df.fragment_id.unique()

    # matrix A is avg fragment
    avg_fragment = avg_fragment[:no_atoms_central]
    A = np.array([np.array(avg_fragment.x), np.array(avg_fragment.y), np.array(avg_fragment.z)]).T
    n = A.shape[0]

    rmse_list = []
    # matrix B is first 100 fragments
    for unique_frag in unique_frags[:100]:
        frag = df[df.fragment_id == unique_frag][:no_atoms_central]
        B = np.array([np.array(frag.x), np.array(frag.y), np.array(frag.z)]).T

        rmse_list.append(calc_rmse(A, B, n))

    avg_rmse = np.mean(rmse_list)
    print(f"Average RMSE average fragment: {avg_rmse :.2f}")

    return avg_rmse


if __name__ == "__main__":
    main()
