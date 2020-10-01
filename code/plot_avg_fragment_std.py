import pandas as pd
import numpy as np
import math

from helpers.headers import COLORS
import matplotlib.pyplot as plt

def main():
    _, ax = plt.subplots()
    plt.title('Mean of fragment centroid\nConfidence Interval (99%)')
    plt.xlabel('Amount of fragments for mean')
    plt.ylabel('Mean')

    ax = plot_ci(color=COLORS[0], ax=ax, filename='results/NO3_CO_vdw.5/NO3_CO_vdw.5_avg_fragment.hdf', label="NO3_CO")
    ax = plot_ci(color=COLORS[1], ax=ax, filename='results/NO3_H2O_vdw.5/NO3_H2O_vdw.5_avg_fragment.hdf', label="NO3_H2O")
    ax = plot_ci(color=COLORS[2], ax=ax, filename='results/NO3_C6H5R_vdw.5/NO3_C6H5R_vdw.5_avg_fragment.hdf', label="NO3_C6H5R")
    # ax = plot_ci(color=COLORS[3], ax=ax, filename='results/ArCI_XH_vdw.5/ArCI_XH_vdw.5_avg_fragment.hdf', label="ArCI_XH")

    plt.legend()
    plt.show()

def plot_ci(color, ax, filename, label):
    df = pd.read_hdf(filename, 'key')

    columns = df.columns
    df['combined'] = df.sum(axis=1)/len(columns)
    # df['combinedy'] = df.sum(axis=1)/len(columns)
    # df['combinedz'] = df.sum(axis=1)/len(columns)

    xdata = np.arange(10, len(df), 10)

    means = []

    # confidence interval of 99%
    ci = []
    z = 2.576

    for xpoint in xdata:
        n = xpoint
        std = float(df[:xpoint]['combined'].std())

        means.append(float(df[:xpoint]['combined'].mean()))

        ci.append((z * std)/math.sqrt(n))
    
    ax.plot(xdata, means, color=color, label=label)
    ax.fill_between(xdata, (np.array(means) - np.array(ci)), (np.array(means) + np.array(ci)), color=color, alpha=0.1)

    return ax


if __name__ == "__main__":
    main()