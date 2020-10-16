import matplotlib.pyplot as plt

from matplotlib.ticker import PercentFormatter


def main():
    with open("dihedrals_dists_angles.csv", 'r') as inputfile:
        lines = inputfile.readlines()

    dists, angles, torsions = [], [], []

    for line in lines[1:]:
        information = line.split(',')

        torsions.append(float(information[2]))
        torsions.append(float(information[3]))
        torsions.append(float(information[4]))

        angles.append(float(information[5]))
        angles.append(float(information[6]))
        angles.append(float(information[7]))

        dists.append(float(information[8]))
        dists.append(float(information[9]))
        dists.append(float(information[10]))

    print('Average bond lenght:', sum(dists)/len(dists))

    plt.figure()
    plt.title("C-H bond lengths in RCOMe")
    plt.hist(dists, bins=200, density=True)
    plt.xlabel("Bond lenght (A)")
    plt.ylabel("Appearance")
    plt.savefig("bondlengths.pdf")

    print('Average C-C-H angle:', sum(angles)/len(angles))

    plt.figure()
    plt.title("C-C-H angles in RCOMe")
    plt.hist(angles, bins=360, density=True)
    plt.xlabel("Angle")
    plt.ylabel("Appearance")
    plt.savefig("angles.pdf")

    plt.figure()
    plt.title("O-C-C-H Dihedrals in RCOMe")
    plt.hist(torsions, bins=360, density=True)
    plt.xlabel(r'Torsion angle ($^{\circ}$)')
    plt.ylabel("Percentage")
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig("torsion.svg")


if __name__ == "__main__":
    main()
