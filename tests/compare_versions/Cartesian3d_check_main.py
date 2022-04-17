import unittest

import GMatElastic.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.hdf5", "r") as data:

            mat = GMat.Array2d(data["/shape"][...])

            iden = data["/model/I"][...]
            K = data["/model/K"][...]
            G = data["/model/G"][...]

            mat.setElastic(iden, K, G)

            for i in range(20):

                GradU = data[f"/random/{i:d}/GradU"][...]

                Eps = np.einsum("...ijkl,...lk->...ij", mat.I4s(), GradU)
                mat.setStrain(Eps)

                self.assertTrue(np.allclose(mat.Stress(), data[f"/random/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.Tangent(), data[f"/random/{i:d}/Tangent"][...]))


if __name__ == "__main__":

    unittest.main()
