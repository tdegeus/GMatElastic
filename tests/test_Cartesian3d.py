import unittest

import GMatElastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import numpy as np


class Test_main(unittest.TestCase):
    """ """

    def test_Epseq_Sigeq(self):

        A = np.zeros((2, 3, 3, 3))
        A[..., 0, 1] = 1
        A[..., 1, 0] = 1

        self.assertTrue(np.allclose(GMat.Epseq(A), 2 / np.sqrt(3) * np.ones(A.shape[:-2])))
        self.assertTrue(np.allclose(GMat.Sigeq(A), np.sqrt(3.0) * np.ones(A.shape[:-2])))

    def test_Elastic(self):

        shape = [2, 3]
        K = np.random.random(shape)
        G = np.random.random(shape)
        mat = GMat.Elastic2d(K, G)

        gamma = np.random.random(shape)
        epsm = np.random.random(shape)

        Eps = np.zeros(shape + [3, 3])
        Eps[..., 0, 0] = epsm
        Eps[..., 1, 1] = epsm
        Eps[..., 2, 2] = epsm
        Eps[..., 0, 1] = gamma
        Eps[..., 1, 0] = gamma
        mat.Eps = Eps

        Sig = np.zeros(shape + [3, 3])
        Sig[..., 0, 0] = 3 * K * epsm
        Sig[..., 1, 1] = 3 * K * epsm
        Sig[..., 2, 2] = 3 * K * epsm
        Sig[..., 0, 1] = 2 * G * gamma
        Sig[..., 1, 0] = 2 * G * gamma

        self.assertTrue(np.allclose(GMat.Epseq(mat.Eps), 2 / np.sqrt(3) * gamma))
        self.assertTrue(np.allclose(GMat.Sigeq(mat.Sig), 2 * np.sqrt(3) * G * gamma))
        self.assertTrue(np.allclose(mat.Sig, Sig))
        self.assertTrue(np.allclose(tensor.A4_ddot_B2(mat.C, mat.Eps), Sig))
        self.assertTrue(np.allclose(mat.energy, 3 * K * epsm**2 + 2 * G * gamma**2))

    def test_tangent(self):

        shape = [2, 3]
        Eps = np.random.random(shape + [3, 3])
        Eps = tensor.A4_ddot_B2(tensor.Array2d(shape).I4s, Eps)
        mat = GMat.Elastic2d(np.random.random(shape), np.random.random(shape))
        mat.Eps = Eps
        self.assertTrue(np.allclose(tensor.A4_ddot_B2(mat.C, mat.Eps), mat.Sig))


if __name__ == "__main__":

    unittest.main()
