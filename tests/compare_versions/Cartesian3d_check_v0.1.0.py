import unittest

import GMatElastic.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.hdf5", "r") as data:

            shape = data["/shape"][...]

            i = np.eye(3)
            I4 = np.einsum("xy,ijkl->xyijkl", np.ones(shape), np.einsum("il,jk", i, i))
            I4rt = np.einsum("xy,ijkl->xyijkl", np.ones(shape), np.einsum("ik,jl", i, i))
            I4s = (I4 + I4rt) / 2.0

            mat = GMat.Matrix(shape[0], shape[1])

            iden = data["/model/I"][...]
            K = data["/model/K"][...]
            G = data["/model/G"][...]

            mat.setElastic(iden, K, G)

            for i in range(20):

                GradU = data[f"/random/{i:d}/GradU"][...]

                Eps = np.einsum("...ijkl,...lk->...ij", I4s, GradU)

                self.assertTrue(np.allclose(mat.Stress(Eps), data[f"/random/{i:d}/Stress"][...]))
                self.assertTrue(
                    np.allclose(mat.Tangent(Eps)[1], data[f"/random/{i:d}/Tangent"][...])
                )


if __name__ == "__main__":

    unittest.main()
