from helpers.headers import CUT_OFF_ZERO
import math
import numpy as np


def compare_distances(oldfragment, newfragment):
    distances1 = calc_distances(oldfragment)
    distances2 = calc_distances(newfragment)

    distances1.sort()
    distances2.sort()

    np.testing.assert_allclose(distances1, distances2, err_msg=oldfragment.from_entry + oldfragment.fragment_id +
                               "after rotating the distances differ")


def calc_distances(fragment):
    distances = []
    counter = 0

    for atom1 in fragment.atoms.values():
        for atom2 in list(fragment.atoms.values())[counter:]:
            if not atom1 == atom2:
                diffx = abs(atom1.x - atom2.x)
                diffy = abs(atom1.y - atom2.y)
                diffz = abs(atom1.z - atom2.z)

                distance = math.sqrt(diffx**2 + diffy**2 + diffz**2)
                distances.append(distance)

        counter += 1

    return distances


def test_centering(fragment, atom_to_center):

    atom = fragment.atoms[atom_to_center]

    assert (atom.x > -CUT_OFF_ZERO and atom.x < CUT_OFF_ZERO and
            atom.y > -CUT_OFF_ZERO and atom.y < CUT_OFF_ZERO and
            atom.z > -CUT_OFF_ZERO and atom.z < CUT_OFF_ZERO), atom_to_center + " atom is not centered right"


def check_new_alignment(fragment, atoms_to_put_in_plane):
    atom1 = atoms_to_put_in_plane[0]

    assert (atom1.x >= 0.0 and atom1.y > -CUT_OFF_ZERO and atom1.y < CUT_OFF_ZERO), fragment.from_entry +\
        " " + atom1.label + " is not in xy plane, check first rotation"

    assert (atom1.x > 0 and atom1.y > -CUT_OFF_ZERO and atom1.y < CUT_OFF_ZERO and atom1.z > -CUT_OFF_ZERO
            and atom1.z < CUT_OFF_ZERO), fragment.from_entry + " " + atom1.label + " is not on the positive x-axis,\
            check second rotation"

    # the second atom is supposed to be on the xy plane, so z has to be zero
    atom2 = atoms_to_put_in_plane[1]

    assert (atom2.z > -CUT_OFF_ZERO and atom2.z < CUT_OFF_ZERO), fragment.from_entry + " " + atom2.label + "\
            is not in the xy plane (check third rotation)"


def test_count(density_df, points_df):
    assert (len(points_df) == (density_df.amount_C.sum() + density_df.amount_O.sum())),\
            "Amount of points is not divided into bins correctly"
