# Implements Kabsch algorithm - best fit.

import numpy as np
from math import sqrt


def rigid_transform_3D(A, B):
    """ Input:
             Nominal  A Nx3 matrix of points
             Measured B Nx3 matrix of points
        Returns s,R,t
            R = 3x3 rotation matrix (B to A)
            t = 3x1 translation vector (B to A)
    """

    assert len(A) == len(B)

    N = A.shape[0]  # total points
    print(N)

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.transpose(BB) * AA

    print(H)

    U, S, Vt = np.linalg.svd(H)
    print("U", U)
    print("S", S)
    print("Vt", Vt)

    rotation_matrix = Vt.T * U.T

    # special reflection case
    if np.linalg.det(rotation_matrix) < 0:
        print("Reflection detected")
        Vt[2, :] *= -1
        rotation_matrix = Vt.T * U.T

    translation_vector = -rotation_matrix * centroid_B.T + centroid_A.T

    return rotation_matrix, translation_vector

A = np.matrix([[ 0.,       0.,       0.     ],
                [-0.74967, -0.49435 , 1.1711 ],
                [-1.71385, -1.22366,  0.97276],
                [-0.37119, -0.16442  ,2.28682]])

B = np.matrix([[ 2.86625,  7.31144,  9.77711],
                [ 1.74101,  7.25375 ,10.72657],
                [ 1.90668 , 6.63789, 11.77807],
                [ 0.701 ,   7.8186 , 10.42211]]) 

B_total = np.matrix([[ 2.86625,  7.31144,  9.77711],
                        [ 1.74101,  7.25375, 10.72657],
                        [ 1.90668,  6.63789, 11.77807],
                        [ 0.701  ,  7.8186 , 10.42211],
                        [ 0.82358, 10.33899, 13.40954],
                        [ 0.29058,  9.76922, 12.79957],
                        [ 2.22175,  9.74507, 13.57359],
                        [ 0.36679, 10.36024 ,14.28762],
                        [ 0.90689, 11.74762, 12.8455 ]])

n = B.shape[0]

# recover the transformation
rotation_matrix, translation_vector = rigid_transform_3D(A, B)

# Find the error
B2 = (rotation_matrix * B.T) + np.tile(translation_vector, (1, n))
B2 = B2.T
err = A - B2
err = np.multiply(err, err)
err = np.sum(err)
rmse = sqrt(err / n)

n = B_total.shape[0]
B2_total = (rotation_matrix * B_total.T) + np.tile(translation_vector, (1, n))

print(n)
print(type(translation_vector))

print(B2_total.T)


print("Rotation")
print(rotation_matrix)

print("Translation")
print(translation_vector)

print("RMSE:", rmse)
print("If RMSE is near zero, the function is correct!")