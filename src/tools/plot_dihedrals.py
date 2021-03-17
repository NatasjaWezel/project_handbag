import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
from matplotlib.ticker import PercentFormatter

df_rcome = pd.read_csv("../data/RCOMe/RCOMe.csv")
df_ret = pd.read_csv("../data/REt/REt.csv")

fig, axs = plt.subplots(2, sharex=True, figsize=(10, 4))
plt.subplots_adjust(bottom=0.15)
plt.suptitle("Torsion angles in REt and RCOMe")


def make_torsionlists(df):
    tor1 = df["TOR1"].tolist()
    tor2 = df["TOR2"].tolist()
    tor3 = df["TOR3"].tolist()

    torsions = tor1 + tor2 + tor3

    return torsions


def average_distance(df):
    dist1 = df["DIST1"].tolist()
    dist2 = df["DIST2"].tolist()
    dist3 = df["DIST3"].tolist()

    distances = dist1 + dist2 + dist3

    return np.mean(distances)


def average_angle(df):
    ang1 = df["ANG1"].tolist()
    ang2 = df["ANG2"].tolist()
    ang3 = df["ANG3"].tolist()

    angles = ang1 + ang2 + ang3

    return np.mean(angles)


print("RCOMe")
print(average_angle(df_rcome))
print(average_distance(df_rcome))

print("REt")
print(average_angle(df_ret))
print(average_distance(df_ret))

torsions_rcome = make_torsionlists(df_rcome)
axs[0].hist(torsions_rcome, bins=360, density=True, color="blue", label="O-C-C-H dihedral RCOMe")
axs[0].grid(True)
axs[0].yaxis.set_major_formatter(PercentFormatter(1))
axs[0].legend(loc='upper right')
axs[0].set_ylabel("Percentage")

torsions_ret = make_torsionlists(df_ret)
axs[1].hist(torsions_ret, bins=360, density=True, color="orange", label="H-C-C-H dihedral REt")
axs[1].grid(True)
axs[1].yaxis.set_major_formatter(PercentFormatter(1))
axs[1].legend(loc='upper right')
axs[1].set_ylabel("Percentage")

plt.xlabel(r'Torsion angle ($^{\circ}$)')

plt.savefig("../../results/torsions.svg", format="svg")
plt.show()
