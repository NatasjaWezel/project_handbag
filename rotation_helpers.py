import numpy as np
import math
import copy

from helpers import plot_fragments

def perform_rotations(fragment, atoms_to_put_in_plane, plot):
    """ Performs three rotations to lie three of the atoms in the xy plane, one of those
        on the x-axis. """ 

    """ First rotation: puts first atom on xy-plane if it already was on a plane, 
        and above the x-axis if it wasn't by rotating around the z-axis. """
    # if first atom doesn't lie in any plane, some extra preparation is required
    atom = fragment.atoms[atoms_to_put_in_plane[0]]
    
    not_in_any_plane = False
    
    if not atom.x == 0.0 and not atom.y == 0.0 and not atom.z == 0.0:
        not_in_any_plane = True
    
    fragment = rotate_fragment(fragment=fragment, atom=atoms_to_put_in_plane[0], ax="z", not_in_any_plane=not_in_any_plane)
   
    """ Second rotation: puts first atom on x-axis by rotating around y-axis. """
    fragment = rotate_fragment(fragment=fragment, atom=atoms_to_put_in_plane[0], ax="y", not_in_any_plane=False)

    """ Third rotation: puts second atom in x-y plane by rotating around the x-axis. """
    fragment = rotate_fragment(fragment=fragment, atom=atoms_to_put_in_plane[1], ax="x", not_in_any_plane=True)

    return fragment


def rotate_fragment(fragment, atom, ax, not_in_any_plane):
    """ Finds coordination vector, calculates its angle with corresponding ax and
        then rotates the molecule. Returns the molecule with it's new coordinates. """ 

    atom = fragment.atoms[atom]
    coord_vector = [atom.x, atom.y, atom.z]

    # if the atom is not in a single plane, project it onto the xy plane for the first rotation
    if not_in_any_plane:
        if ax == "x":
            coord_vector = [0, atom.y, atom.z]
        elif ax == "y":
            coord_vector = [atom.x, atom.y, atom.z]
        else:
            coord_vector = [atom.x, atom.y, 0]
        
    # calculate angles and perform rotation
    alpha, beta, _ = find_angles(coord_vector=coord_vector)
    
    if ax == "x":
        angle = beta

        if atom.z < 0:
            angle = -beta
    elif ax == "y":
        angle = -alpha

        if atom.z < 0:
            alpha = alpha

    elif ax == "z":
        angle = alpha

        if atom.y < 0:
            angle = -alpha

    fragment = calculate_rotation(fragment=fragment, angle=angle, ax=ax)

    return fragment

def find_angles(coord_vector):
    """ Rotates the molecule so that the contact fragment is always in the same position. 
        Rotates only the important part of the molecule. """

    # x, y, z unitary row vectors:
    x, y, z =  np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])

    point_vector = np.array(coord_vector)
    
    # formula: u.v = |u|.|v|.cos(alpha)
    # alpha = arccos((u.v)/(|u|.|v|))
    alpha = np.arccos(np.dot(point_vector, x) / np.linalg.norm(point_vector))
    beta = np.arccos(np.dot(point_vector, y) / np.linalg.norm(point_vector))
    gamma = np.arccos(np.dot(point_vector, z) / np.linalg.norm(point_vector))

    return alpha, beta, gamma


def calculate_rotation(fragment, angle, ax):
    # rotate all coordinates according to the previously defined rotation matrices
    for atom in fragment.atoms.values():
        coord_vector = np.array([atom.x, atom.y, atom.z])

        if ax == "x":
            coord_vector = np.dot(coord_vector, rotate_x(angle))
        elif ax == "y":
            coord_vector = np.dot(coord_vector, rotate_y(angle))
        else:  
            coord_vector = np.dot(coord_vector, rotate_z(angle))
            
        atom.x, atom.y, atom.z = coord_vector[0], coord_vector[1], coord_vector[2]

    return fragment


def rotate_x(angle):
    """ Rotation matrix for rotation around x-axis. """

    return np.array(( [1,             0,              0], 
                        [0,             np.cos(angle),  -np.sin(angle)], 
                        [0,             np.sin(angle),  np.cos(angle)]))

def rotate_y(angle):
    """ Rotation matrix for rotation around y-axis. """

    return np.array(( [np.cos(angle),  0,              -np.sin(angle)], 
                        [0,             1,              0], 
                        [np.sin(angle),  0,              np.cos(angle)]))

def rotate_z(angle):
    """ Rotation matrix for rotation around z-axis. """

    return np.array(( [np.cos(angle), -np.sin(angle), 0], 
                        [np.sin(angle), np.cos(angle),  0], 
                        [0,             0,              1]))