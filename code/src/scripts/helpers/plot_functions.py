import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D


def plot_density(ax, df, settings, lower_lim=0.0001):
    df['ymiddle'] = (df['ystart'] * 2 + settings.resolution) / 2
    df['xmiddle'] = (df['xstart'] * 2 + settings.resolution) / 2
    df['zmiddle'] = (df['zstart'] * 2 + settings.resolution) / 2

    # normalize column
    df[settings.to_count_contact + "_normalized"] = df[settings.to_count_contact] / df[settings.to_count_contact].sum()

    points = df[df[settings.to_count_contact + "_normalized"] > lower_lim]

    norm = plt.Normalize(lower_lim, points[settings.to_count_contact + "_normalized"].max())
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue", "fuchsia", "red"])

    p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle),
                   s=list(10000 * points[settings.to_count_contact + "_normalized"]),
                   c=list(points[settings.to_count_contact + "_normalized"]),
                   cmap=cmap,
                   norm=norm)

    ax.set_xlabel('X axis ($\\AA$)')
    ax.set_ylabel('Y axis ($\\AA$)')
    ax.set_zlabel('Z axis ($\\AA$)')

    return p, ax


def plot_fragment_colored(ax, fragment):
    for _, atom in fragment.iterrows():
        color = get_atom_color(atom)

        ax.scatter(atom.x, atom.y, atom.z, c=color, s=atom.cov_radius*500, edgecolor="black")


def get_atom_color(atom_row):
    if atom_row.symbol == "O":
        return "red"
    elif atom_row.symbol == "H":
        return "white"
    elif atom_row.symbol == "N":
        return "blue"
    elif atom_row.symbol == "C":
        return "black"
    elif atom_row.symbol == "I":
        return "yellow"
    elif "aH" in atom_row.label:
        return "Fuchsia"
    else:
        return "Pink"


def plot_fragments(df, amount, COLORS):
    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    i = 0
    unique = list(df.fragment_id.unique())

    for _id in unique[:amount]:
        fragment = df[df.fragment_id == _id]

        color = COLORS[i % len(COLORS)]
        ax.scatter(fragment.x, fragment.y, fragment.z, color=color, label=_id)

        i += 1

    ax.legend()
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


def plot_vdw_spheres(avg_fragment, ax, extra=0):
    spheres = []

    for _, atom in avg_fragment.iterrows():
        r = atom.vdw_radius + extra

        theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

        x = r * np.sin(phi) * np.cos(theta) + atom.x
        y = r * np.sin(phi) * np.sin(theta) + atom.y
        z = r * np.cos(phi) + atom.z

        color = get_atom_color(atom)
        sphere = ax.plot_surface(x, y, z, color=color, alpha=0.5, linewidth=0)
        spheres.append(sphere)

    assert len(spheres) == len(avg_fragment), "Something went wrong with plotting the vdw surfaces"

    return ax, spheres
