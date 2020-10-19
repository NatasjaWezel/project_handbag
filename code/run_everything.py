###
import os

import numpy as np


def main():
    # central_groups = ["ArCI", "H2O", "NO3", "RC6F6", "RC6H5", "RCOMe", "RNO2"]
    central_groups = ["H2O"]
    contact_groups = ["ArCH", "CF", "R2CO", "RC6H5", "RCN", "XH"]
    resolutions = np.arange(0.2, 1.6, 0.1)

    for central_group in central_groups[:1]:
        for contact_group in contact_groups:
            datafile = ".\\data\\" + central_group + "\\" + central_group + "_" + contact_group + "_vdw.5.cor"
            result1 = ".\\results\\" + central_group + "\\" + central_group + "_" + contact_group + "_vdw.5\\"\
                      + central_group + "_" + contact_group + "_vdw.5_aligned.csv"

            os.system("python 1_load_from_coords.py " + datafile)

            for resolution in resolutions:
                print("\nCalculating density for central group: ", central_group, " contact group: ", contact_group,
                      "resolution: ", str(round(resolution, 2)))

                os.system("python 4_calc_density.py " + result1 + " " + str(round(resolution, 2)) + " " + "centroid")

                os.system("python 5_analyze_density.py " + result1 + " " + str(round(resolution, 2)) + " " + "centroid")


if __name__ == "__main__":
    main()
