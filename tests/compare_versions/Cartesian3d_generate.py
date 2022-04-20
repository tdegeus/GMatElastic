import GMatElastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import h5py
import numpy as np

with h5py.File("Cartesian3d_random.h5", "w") as data:

    shape = [1000, 4]
    K = 0.5 + np.random.random(shape)
    G = 0.5 + np.random.random(shape)

    data["/model/K"] = K
    data["/model/G"] = G

    mat = GMat.Elastic2d(K, G)
    I4s = tensor.Array2d(mat.shape).I4s

    for i in range(20):

        GradU = 200 * np.random.random(shape + [3, 3])
        data[f"/random/{i:d}/GradU"] = GradU
        mat.Eps = tensor.A4_ddot_B2(I4s, GradU)

        data[f"/random/{i:d}/Stress"] = mat.Sig
        data[f"/random/{i:d}/Tangent"] = mat.C
