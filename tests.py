from headers import *

def test_centering(fragment, atom_to_center):

    atom = fragment.atoms[atom_to_center]

    assert (atom.x > -CUT_OFF_ZERO and atom.x < CUT_OFF_ZERO and 
                atom.y > -CUT_OFF_ZERO and atom.y < CUT_OFF_ZERO and
                atom.z > -CUT_OFF_ZERO and atom.z < CUT_OFF_ZERO), atom_to_center + " atom is not centered right"

def test_rotation1(fragment, atoms_to_put_in_plane):
    atom1 = fragment.atoms[atoms_to_put_in_plane[0]]
    
    assert (atom1.x >= 0.0 and atom1.y > -CUT_OFF_ZERO and atom1.y < CUT_OFF_ZERO), atom1.label + " is not in xy plane, check first rotation"

def test_rotation2(fragment, atoms_to_put_in_plane):
    # the first atom is supposed to be on the x axis, so y and z have to be 0
    atom1 = fragment.atoms[atoms_to_put_in_plane[0]]
    
    assert (atom1.x > 0 and atom1.y > -CUT_OFF_ZERO and atom1.y < CUT_OFF_ZERO and atom1.z > -CUT_OFF_ZERO and atom1.z < CUT_OFF_ZERO), atom1.label + " is not on the positive x-axis, check second rotation"

def test_rotation3(fragment, atoms_to_put_in_plane):
    # the second atom is supposed to be on the xy plane, so z has to be zero
    atom2 = fragment.atoms[atoms_to_put_in_plane[1]]

    assert (atom2.z > -CUT_OFF_ZERO and atom2.z < CUT_OFF_ZERO), atom2.label + " is not in the xy plane (check third rotation)"

    