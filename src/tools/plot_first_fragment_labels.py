import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.append('../scripts')

from classes.Radii import Radii
from classes.Settings import AlignmentSettings


def main():
    datafile = sys.argv[1]

    settings = AlignmentSettings('../', datafile)

    plot_fragment_with_labels(settings)


def plot_fragment_with_labels(settings):

    fp = open(settings.label_data)
    labels = fp.readline().strip().split(',')
    atoms = fp.readline().strip().split(',')
    fp.close()

    to_delete = []
    for i in range(len(labels)):
        if "LAB" not in labels[i]:
            to_delete.append(i)

    for i in reversed(to_delete):
        del labels[i]
        del atoms[i]

    fp = open(settings.coordinate_file)
    line = fp.readline()
    line = fp.readline()

    firstfragment = True
    dictionary = {}

    while firstfragment:
        information = line.split()
        x, y, z = float(information[1]), float(information[2]), float(information[3])

        if information[0] in atoms:
            dictionary[information[0]] = [x, y, z]

        line = fp.readline()

        if "**FRAG**" in line:
            firstfragment = False

    xs, ys, zs = [], [], []
    for key, value in dictionary.items():
        xs.append(value[0])
        ys.append(value[1])
        zs.append(value[2])

    centroid = [np.mean(xs), np.mean(ys), np.mean(zs)]

    for key, value in dictionary.items():
        dictionary[key][0] = dictionary[key][0] - centroid[0]
        dictionary[key][1] = dictionary[key][1] - centroid[1]
        dictionary[key][2] = dictionary[key][2] - centroid[2]

    fp.close()

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    for i, atom in enumerate(atoms):
        x, y, z = dictionary[atom][0], dictionary[atom][1], dictionary[atom][2]

        if 'H' in atom:
            color = 'grey'
            label = 'H'
        elif 'O' in atom:
            color = 'red'
            label = 'O'
        elif 'F' in atom:
            color = 'orchid'
            label = 'F'
        elif 'N' in atom:
            color = 'blue'
            label = 'N'
        elif 'C' in atom:
            color = 'black'
            label = 'C'
        else:
            color = "purple"
            label = "other then H, O, F, N or C"

        ax.scatter(x, y, z, s=100, edgecolors="black", color=color, label=label)
        ax.text(x+0.1, y+0.1, z+0.1, labels[i])

    for i, atom in enumerate(atoms):
        x1, y1, z1 = dictionary[atom][0], dictionary[atom][1], dictionary[atom][2]
        for i2, atom2 in enumerate(atoms[i:]):
            x2, y2, z2 = dictionary[atom2][0], dictionary[atom2][1], dictionary[atom2][2]
            distance = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5

            radii = Radii('../files/radii.csv')

            if distance < (radii.get_cov_radius(atom[:1]) + radii.get_cov_radius(atom2[:1]) + 0.01):
                plt.plot([x1, x2], [y1, y2], [z1, z2], color='grey')

    # make cubic
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()

    lim_min, lim_max = min(xlim[0], ylim[0], zlim[0]), max(xlim[1], ylim[1], zlim[1])

    ax.set_xlim((lim_min, lim_max))
    ax.set_ylim((lim_min, lim_max))
    ax.set_zlim((lim_min, lim_max))

    # Hide grid lines
#     ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

#     ax.axis('off')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.title(settings.central_name + "-" + settings.contact_name)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

    plt.show()


if __name__ == "__main__":
    main()
