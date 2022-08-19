import unittest

import GMatElastic.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.h5") as data:

            mat = GMat.Elastic2d(data["/K"][...], data["/G"][...])

            for i in range(20):

                mat.Eps = data[f"/data/{i:d}/Eps"][...]
                self.assertTrue(np.allclose(mat.Sig, data[f"/data/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.C, data[f"/data/{i:d}/Tangent"][...]))


if __name__ == "__main__":

    unittest.main()
