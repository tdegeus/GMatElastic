import h5py
import numpy as np
import GMatElastic.Cartesian3d as GMat

with h5py.File('Cartesian3d_random.hdf5', 'w') as data:

    nelem = 1000
    nip = 4

    shape = np.array([nelem, nip], np.int)

    data['/shape'] = shape

    mat = GMat.Array2d(shape)

    I = np.ones(shape)
    K = np.ones(n)
    G = np.ones(n)

    data['/elastic/I'] = I
    data['/elastic/K'] = K
    data['/elastic/G'] = G

    mat.setElastic(I, K, G)
    mat.check()

    for i in range(20):

        GradU = 200 * np.random.random([nelem, nip, 2, 2])

        data['/random/{0:d}/GradU'.format(i)] = GradU

        Eps = np.einsum('...ijkl,...lk->...ij', mat.I4s(), GradU)
        mat.setStrain(Eps)

        data['/random/{0:d}/Stress'.format(i)] = mat.Stress()
        data['/random/{0:d}/Tangent'.format(i)] = mat.Tangent()
