import unittest

import GMatElastic.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.h5") as data:

            K = data["/model/K"][...]
            G = data["/model/G"][...]
            mat = GMat.Array2d(K.shape)

            for i in range(K.shape[0]):
                for j in range(K.shape[1]):
                    iden = np.zeros(K.shape, bool)
                    iden[i, j] = True
                    mat.setElastic(iden, K[i, j], G[i, j])

            for i in range(20):

                GradU = data[f"/random/{i:d}/GradU"][...]

                Eps = np.einsum("...ijkl,...lk->...ij", mat.I4s(), GradU)
                mat.setStrain(Eps)

                self.assertTrue(np.allclose(mat.Stress(), data[f"/random/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.Tangent(), data[f"/random/{i:d}/Tangent"][...]))


if __name__ == "__main__":

    unittest.main()
