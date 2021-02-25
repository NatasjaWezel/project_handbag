import sys

import pandas as pd
import numpy as np

# example for input avg fragment rcome
ANGLES = ['LAB2-LAB3-LAB5', 'LAB2-LAB3-LAB6', 'LAB2-LAB3-LAB7', 'LAB2-LAB3-aH1', 'LAB2-LAB3-aH2']
DISTANCES = ['LAB3-LAB5', 'LAB3-LAB6', 'LAB3-LAB7', 'LAB3-aH1', 'LAB3-aH2']

# example for input avg fragment ret
# ANGLES = ['LAB5-LAB1-LAB4', 'LAB5-LAB1-LAB2', 'LAB5-LAB1-LAB3', 'LAB5-LAB1-aH1', 'LAB5-LAB1-aH2']
# DISTANCES = ['LAB1-LAB2', 'LAB1-LAB3', 'LAB1-LAB4', 'LAB1-aH1', 'LAB1-aH2']


def main():

    if len(sys.argv) != 2:
        print("Usage: python calc_angles.py <path/to/inputfile>")
        print("Don't forget to specify the atoms")
        sys.exit(1)

    df = pd.read_csv(sys.argv[1])
    print(df)

    for angle in ANGLES:
        atom1, atom2, atom3 = angle.split('-')
        v1 = np.array((df[df.label == atom1].x, df[df.label == atom1].y, df[df.label == atom1].z))
        v2 = np.array((df[df.label == atom2].x, df[df.label == atom2].y, df[df.label == atom2].z))
        v3 = np.array((df[df.label == atom3].x, df[df.label == atom3].y, df[df.label == atom3].z))

        # print(v1, v2, v2-v1)

        degrees = calc_angle(v2-v1, v3-v2)
        print(f"Angle {angle}: {float(degrees) :.2f}")

    for distance in DISTANCES:
        atom1, atom2 = distance.split('-')
        v1 = np.array((df[df.label == atom1].x, df[df.label == atom1].y, df[df.label == atom1].z))
        v2 = np.array((df[df.label == atom2].x, df[df.label == atom2].y, df[df.label == atom2].z))

        length = calc_distance(v1, v2)
        print(f"Distance {distance}: {float(length) :.2f}")



def calc_distance(v1, v2):
    return np.linalg.norm(v1 - v2)


def calc_angle(v1, v2):
    return 180 - np.rad2deg(np.arccos(np.dot(v1.T, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))


if __name__ == "__main__":
    main()
