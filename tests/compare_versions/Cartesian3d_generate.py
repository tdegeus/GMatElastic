import GMatElastic.Cartesian3d as GMat
import h5py
import numpy as np

with h5py.File("Cartesian3d_random.hdf5", "w") as data:

    nelem = 1000
    nip = 4

    shape = np.array([nelem, nip], np.int)

    data["/shape"] = shape

    mat = GMat.Array2d(shape)

    iden = np.ones(shape).astype(np.int)
    n = iden.size
    K = 12.3
    G = 45.6

    data["/model/I"] = iden
    data["/model/K"] = K
    data["/model/G"] = G

    mat.setElastic(iden, K, G)

    for i in range(20):

        GradU = 200 * np.random.random([nelem, nip, 3, 3])

        data[f"/random/{i:d}/GradU"] = GradU

        Eps = np.einsum("...ijkl,...lk->...ij", mat.I4s(), GradU)
        mat.setStrain(Eps)

        data[f"/random/{i:d}/Stress"] = mat.Stress()
        data[f"/random/{i:d}/Tangent"] = mat.Tangent()
