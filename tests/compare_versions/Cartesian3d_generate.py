import GMatElastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import h5py
import numpy as np

with h5py.File("Cartesian3d_random.h5", "w") as data:

    shape = [1000, 4]

    mat = GMat.Elastic2d(
        K=0.5 + np.random.random(shape),
        G=0.5 + np.random.random(shape),
    )

    data["K"] = mat.K
    data["G"] = mat.G

    for i in range(20):

        mat.Eps = tensor.Sym(20 * np.random.random(shape + [3, 3]))

        data[f"/data/{i:d}/Eps"] = mat.Eps
        data[f"/data/{i:d}/Stress"] = mat.Sig
        data[f"/data/{i:d}/Tangent"] = mat.C
