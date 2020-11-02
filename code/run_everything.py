###
import os

import numpy as np


def main():
    central_groups = ["NO3"]    # "ArCI", "H2O", "NO3", "RC6F6", "RC6H5", "RCOMe", "RNO2"
    contact_groups = ["XH"]     # "CF", "R2CO", "RC6H5", "RCN", "XH"
    to_count = ["H"]            # , "O", "centroid", "N", "H"
    resolutions = np.arange(0.1, 1.1, 0.1)

    for central_group in central_groups:
        for to_count_contact, contact_group in zip(to_count, contact_groups):
            datafile = ".\\data\\" + central_group + "\\" + central_group + "_" + contact_group + "_vdw.5.cor"
            result1 = ".\\results\\" + central_group + "\\" + central_group + "_" + contact_group + "_vdw.5\\"\
                      + central_group + "_" + contact_group + "_vdw.5_aligned.csv"

            print(datafile)
            print(result1)
            os.system("python 1_load_from_coords.py " + datafile)

            for resolution in resolutions:
                print("\nCalculating density for central group: ", central_group, " contact group: ", contact_group,
                      "resolution: ", str(round(resolution, 2)))

                os.system("python 4_calc_density.py " + result1 + " " + str(round(resolution, 2)) +
                          " " + to_count_contact)

                # os.system("python analyze_density.py " + result1 + " " + str(round(resolution, 2)) +
                #           " " + to_count_contact)


if __name__ == "__main__":
    main()
