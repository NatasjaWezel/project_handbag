import numpy as np

def rotate_molecule(molecule, atom, ax, not_in_any_plane):
    """ Finds coordination vector, calculates its angle with corresponding ax and
        then rotates the molecule. Returns the molecule with it's new coordinates. """ 

    coord_vector = [molecule.highlighted_atoms[atom].x, molecule.highlighted_atoms[atom].y, molecule.highlighted_atoms[atom].z]

    # if the atom is not in a single plane, project it onto the xy plane for the first rotation
    if not_in_any_plane:
        if ax == "x":
            pass
        elif ax == "y":
            pass
        else:
            coord_vector = [molecule.highlighted_atoms[atom].x, molecule.highlighted_atoms[atom].y, 0]
        
    # calculate angles and perform rotation
    alpha, beta, _ = find_angles(coord_vector=coord_vector)
    
    if ax == "x":
        angle = beta
    elif ax == "y":
        angle = -alpha

        if molecule.highlighted_atoms[atom].z < 0:
            alpha = alpha
    else:
        angle = alpha

        if molecule.highlighted_atoms[atom].y < 0:
            angle = -alpha
    
    molecule = calculate_rotation(molecule=molecule, angle=angle, ax=ax)

    return molecule

def find_angles(coord_vector):
    """ Rotates the molecule so that the contact group is always in the same position. 
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


def calculate_rotation(molecule, angle, ax):
    # rotate all coordinates according to the previously defined rotation matrices
    for atom in molecule.highlighted_atoms.values():
        coord_vector = np.array([atom.x, atom.y, atom.z])

        if ax == "x":
            coord_vector = np.dot(coord_vector, rotate_x(angle))
        elif ax == "y":
            coord_vector = np.dot(coord_vector, rotate_y(angle))
        else:  
            coord_vector = np.dot(coord_vector, rotate_z(angle))
            
        atom.x, atom.y, atom.z = coord_vector[0], coord_vector[1], coord_vector[2]

    return molecule


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