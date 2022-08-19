import unittest

import GMatElastic.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.h5") as data:

            K = data["K"][...]
            G = data["G"][...]
            mat = GMat.Array2d(K.shape)

            for i in range(K.shape[0]):
                for j in range(K.shape[1]):
                    iden = np.zeros(K.shape, bool)
                    iden[i, j] = True
                    mat.setElastic(iden, K[i, j], G[i, j])

            for i in range(20):

                Eps = data[f"/data/{i:d}/Eps"][...]
                mat.setStrain(Eps)

                self.assertTrue(np.allclose(mat.Stress(), data[f"/data/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.Tangent(), data[f"/data/{i:d}/Tangent"][...]))


if __name__ == "__main__":

    unittest.main()
