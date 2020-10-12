import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib
from helpers.headers import COLORS

import pandas as pd
import numpy as np

def plot_density(ax, df, resolution):
    df['ymiddle'] = (df['ystart'] * 2 + resolution) / 2
    df['xmiddle'] = (df['xstart'] * 2 + resolution) / 2
    df['zmiddle'] = (df['zstart'] * 2 + resolution) / 2

    # for now only use first column
    column_name = [i for i in df.columns if "amount" in i][0]
    # normalize per column 
    df[column_name] = df[column_name] / df[column_name].sum()

    points = df[df[column_name] > 0.00001]

    norm = plt.Normalize(0.00001, points[column_name].max())
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightblue","fuchsia","red"])
    
    p = ax.scatter(list(points.xmiddle), list(points.ymiddle), list(points.zmiddle), s=list(10000 * points[column_name]), c=list(points[column_name]), cmap=cmap, norm=norm)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    return p, ax


def plot_fragment_colored(ax, fragment):
    if isinstance(fragment, pd.DataFrame):
        return plot_fragment_from_df(ax, fragment)
    else:
        return plot_fragment(ax, fragment)

def plot_fragment_from_df(ax, fragment_df):
    for _, atom in fragment_df.iterrows():
        if "O" in atom.atom_label:
            ax.scatter(atom.atom_x, atom.atom_y, atom.atom_z, c="red", s=30, edgecolor="black")
        elif "N" in atom.atom_label:
            ax.scatter(atom.atom_x, atom.atom_y, atom.atom_z, c="blue", s=30, edgecolor="black")
        elif "C" in atom.atom_label:
            ax.scatter(atom.atom_x, atom.atom_y, atom.atom_z, c="black", s=30, edgecolor="black")
        elif "I" in atom.atom_label:
            ax.scatter(atom.atom_x, atom.atom_y, atom.atom_z, c="orchid", s=30, edgecolor="black")
        else:
            ax.scatter(atom.atom_x, atom.atom_y, atom.atom_z, c="pink", s=30, edgecolor="black")

    return ax        

def plot_fragment(ax, fragment):
    # plot the (average of the) central group 
    for atom in fragment.atoms.values():
        if "O" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="red", s=30, edgecolor="black")
        elif "N" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="blue", s=30, edgecolor="black")
        elif "C" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="black", s=30, edgecolor="black")
        elif "I" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="orchid", s=30, edgecolor="black")
        else:
            ax.scatter(atom.x, atom.y, atom.z, c="pink", s=30, edgecolor="black")

    return ax


def plot_fragments(df, amount):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    i = 0
    unique = list(df.id.unique())

    for _id in unique[:amount]:
        fragment = df[df.id == _id]

        color = COLORS[i % len(COLORS)]
        ax.scatter(fragment.atom_x, fragment.atom_y, fragment.atom_z, color=color, label=_id)

        i += 1

    ax.legend()
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    
    plt.show()

def plot_vdw_spheres(avg_fragment, ax, color, extra=0):
    spheres = []

    for _, atom in avg_fragment.iterrows():
        r = atom.vdw_radius + extra

        theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

        x = r * np.sin(phi) * np.cos(theta) + atom.atom_x
        y = r * np.sin(phi) * np.sin(theta) + atom.atom_y
        z = r * np.cos(phi) + atom.atom_z

        sphere = ax.plot_surface(x, y, z, color=color, alpha=0.2, linewidth=0)
        spheres.append(sphere)

    assert len(spheres) == len(avg_fragment), "Something went wrong with plotting the vdw surfaces"

    return ax, spheres