import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib

def plot_density(ax, to_count, df, resolution):
    df['ymiddle'] = (df['ystart'] + df['yend']) / 2
    df['xmiddle'] = (df['xstart'] + df['xend']) / 2
    df['zmiddle'] = (df['zstart'] + df['zend']) / 2

    # for now only use first column
    column = to_count[0]
    columname = "amount_" + column
    # normalize per column 
    df[columname] = df[columname] / df[columname].sum()

    points = df[df[columname] > 0.001]

    norm = plt.Normalize(0.001, points[columname].max())
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue","fuchsia","red"])
    
    # TODO: fix sizes of points
    p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle), s=list(10000 * points[columname]), c=list(points[columname]), cmap=cmap, norm=norm)

    ax.set_title("4D density plot\n Resolution: " + str(resolution))

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    return p, ax


def plot_fragment_colored(ax, fragment):
    # plot the (average of the) central group 
    for atom in fragment.atoms.values():
        if "O" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="red", s=20, edgecolor="black")
        elif "N" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="blue", s=20, edgecolor="black")

    return ax


def plot_fragments(fragments):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    colors = ["red", "sandybrown", "gold", "chartreuse", "green", 
                "mediumturquoise", "dodgerblue", "darkblue", "slateblue",
                "mediumorchid", "fuchsia"]

    for i, fragment in enumerate(fragments):
        fragment.color = colors[i % len(colors)]

        # plot first atom of the fragment and label it
        atom = list(fragment.atoms.values())[0]
        ax.scatter(atom.x, atom.y, atom.z, color=fragment.color, label=fragment.from_entry + str(fragment.fragment_id))
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black') 
        
        ax = plot_atoms_bonds(ax=ax, fragment=fragment)

    ax.legend()
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    
    plt.show()

def plot_atoms_bonds(ax, fragment):
    """ This function plots the atoms and bonds of a fragment. """ 

    for atom in list(fragment.atoms.values())[1:]:
        ax.scatter(atom.x, atom.y, atom.z, color=fragment.color)
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black')                 
    
    for bond in fragment.bonds:
        # TODO: look into this cuz fragments.atoms.keys are supposed to be labels and bond[0] is supposed to be an atom
        if bond[0] in fragment.atoms.keys() and bond[1] in fragment.atoms.keys():
            x = [fragment.atoms[bond[0]].x, fragment.atoms[bond[1]].x]
            y = [fragment.atoms[bond[0]].y, fragment.atoms[bond[1]].y]
            z = [fragment.atoms[bond[0]].z, fragment.atoms[bond[1]].z]

            ax.plot(x, y, z, color=fragment.color)
    
    return ax