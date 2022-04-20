import unittest

import GMatElastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.h5") as data:

            mat = GMat.Elastic2d(data["/model/K"][...], data["/model/G"][...])
            I4s = tensor.Array2d(mat.shape).I4s

            for i in range(20):

                GradU = data[f"/random/{i:d}/GradU"][...]
                mat.Eps = tensor.A4_ddot_B2(I4s, GradU)
                self.assertTrue(np.allclose(mat.Sig, data[f"/random/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.C, data[f"/random/{i:d}/Tangent"][...]))


if __name__ == "__main__":

    unittest.main()
