from helpers.headers import CUT_OFF_ZERO

def test_centering(fragment, atom_to_center):

    atom = fragment.atoms[atom_to_center]

    assert (atom.x > -CUT_OFF_ZERO and atom.x < CUT_OFF_ZERO and 
                atom.y > -CUT_OFF_ZERO and atom.y < CUT_OFF_ZERO and
                atom.z > -CUT_OFF_ZERO and atom.z < CUT_OFF_ZERO), atom_to_center + " atom is not centered right"

def test_rotation1(fragment, atoms_to_put_in_plane):
    atom1 = atoms_to_put_in_plane[0]
    
    assert (atom1.x >= 0.0 and atom1.y > -CUT_OFF_ZERO and atom1.y < CUT_OFF_ZERO), fragment.from_entry +  " "  + atom1.label + " is not in xy plane, check first rotation"

def test_rotation2(fragment, atoms_to_put_in_plane):
    # the first atom is supposed to be on the x axis, so y and z have to be 0
    atom1 = atoms_to_put_in_plane[0]
    
    assert (atom1.x > 0 and atom1.y > -CUT_OFF_ZERO and atom1.y < CUT_OFF_ZERO and atom1.z > -CUT_OFF_ZERO and atom1.z < CUT_OFF_ZERO), fragment.from_entry +  " "  + atom1.label + " is not on the positive x-axis, check second rotation"

def test_rotation3(fragment, atoms_to_put_in_plane):
    # the second atom is supposed to be on the xy plane, so z has to be zero
    atom2 = atoms_to_put_in_plane[1]

    assert (atom2.z > -CUT_OFF_ZERO and atom2.z < CUT_OFF_ZERO), fragment.from_entry + " "  + atom2.label + " is not in the xy plane (check third rotation)"


def test_count(density_df, points_df):
    assert (len(points_df) == (density_df.amount_C.sum() + density_df.amount_O.sum())), "Amount of points is not divided into bins correctly"
