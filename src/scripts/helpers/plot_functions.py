import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D


def plot_density(ax, df, settings):
    df['ymiddle'] = (df['ystart'] * 2 + settings.resolution) / 2
    df['xmiddle'] = (df['xstart'] * 2 + settings.resolution) / 2
    df['zmiddle'] = (df['zstart'] * 2 + settings.resolution) / 2

    # normalize column
    df[settings.contact_rp + "_normalized"] = df[settings.contact_rp] / df[settings.contact_rp].sum()

    # use threshold to determine lower limit
    lower_lim = settings.threshold * df[settings.contact_rp + "_normalized"].max()
    points = df[df[settings.contact_rp + "_normalized"] > lower_lim]

    norm = plt.Normalize(lower_lim, points[settings.contact_rp + "_normalized"].max())
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue", "fuchsia", "red"])

    p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle),
                   s=list(10000 * points[settings.contact_rp + "_normalized"]),
                   c=list(points[settings.contact_rp + "_normalized"]),
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

    return ax


def get_atom_color(atom_row):
    if atom_row.symbol == "O":
        return "red"
    elif "aH" in atom_row.label:
        return "Fuchsia"
    elif atom_row.symbol == "H":
        return "white"
    elif atom_row.symbol == "N":
        return "blue"
    elif atom_row.symbol == "C":
        return "black"
    elif atom_row.symbol == "I":
        return "yellow"
    else:
        return "Pink"


def plot_fragments(df, amount, COLORS):
    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    plt.title("Superimposed fragments")
    i = 0
    unique = list(df.fragment_id.unique())

    for _id in unique[:amount]:
        fragment = df[df.fragment_id == _id]

        color = COLORS[i % len(COLORS)]

        central = fragment[fragment.label != "-"]
        contact = fragment[fragment.label == "-"]

        ax.scatter(central.x, central.y, central.z, color=color, s=100, label=_id)
        if len(contact) > 0:
            ax.scatter(contact.x, contact.y, contact.z, color=color)

        i += 1

    ax.legend()

    xlim, ylim, zlim = list(ax.get_xlim()), list(ax.get_ylim()), list(ax.get_zlim())
    minn = min([xlim[0], ylim[0], zlim[0]])
    maxx = max([xlim[1], ylim[1], zlim[1]])

    ax.set_xlim((minn, maxx))
    ax.set_ylim((minn, maxx))
    ax.set_zlim((minn, maxx))

    ax.set_xlabel('X axis ($\\AA$)')
    ax.set_ylabel('Y axis ($\\AA$)')
    ax.set_zlabel('Z axis ($\\AA$)')

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
