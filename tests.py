from headers import *

def check_new_fragment_alignment(fragment, atom_to_center, atoms_to_put_in_plane):
    
    test_centering(fragment, atom_to_center)
    test_axis_alignment(fragment, atoms_to_put_in_plane)
    test_plane_alignment(fragment, atoms_to_put_in_plane)

    print("Passed all checks, fragment's new alignment is correct.")

def test_centering(fragment, atom_to_center):
    
    atom = fragment.atoms[atom_to_center]
    if atom.x != 0.0 or atom.y != 0.0 or atom.z != 0.0:
        assert atom_to_center + "atom is not centered right"

def test_axis_alignment(fragment, atoms_to_put_in_plane):
    # the first atom is supposed to be on the x axis, so y and z have to be 0
    atom1 = fragment.atoms[atoms_to_put_in_plane[0]]
    if atom1.y < -CUT_OFF_ZERO or atom1.y > CUT_OFF_ZERO or atom1.z < -CUT_OFF_ZERO or atom1.z > CUT_OFF_ZERO:
        assert atom1.label + " is not on the x axis, check first and second rotations"

def test_plane_alignment(fragment, atoms_to_put_in_plane):
    # the second atom is supposed to be on the xy plane, so z has to be zero
    atom2 = fragment.atoms[atoms_to_put_in_plane[1]]
    if atom2.z < -CUT_OFF_ZERO or atom2.z > CUT_OFF_ZERO:
        assert atom2.label + " is not in the xy plane (check third rotation)"

    