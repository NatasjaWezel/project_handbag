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
    
    # TODO: fix sizes of points
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


def plot_fragments(fragments):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for i, fragment in enumerate(fragments):
        fragment.color = COLORS[i % len(COLORS)]

        atomxs = [atom.x for atom in fragment.atoms.values()]
        atomys = [atom.y for atom in fragment.atoms.values()]
        atomzs = [atom.z for atom in fragment.atoms.values()]

        ax.scatter(atomxs, atomys, atomzs, color=fragment.color, label=fragment.id)

        for atom in fragment.atoms.values():
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


def plot_vdw_spheres(avg_fragment, ax):
    spheres = []

    for atom in avg_fragment.atoms.values():
        r = atom.vdw_radius

        theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

        x = r * np.sin(phi) * np.cos(theta) + atom.x
        y = r * np.sin(phi) * np.sin(theta) + atom.y
        z = r * np.cos(phi) + atom.z

        sphere = ax.plot_surface(x, y, z, color='pink', alpha=0.2, linewidth=0)
        spheres.append(sphere)

    assert len(spheres) == len(avg_fragment.atoms.keys()), "Something went wrong with plotting the vdw surfaces"

    return ax, spheres