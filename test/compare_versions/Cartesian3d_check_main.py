import h5py
import numpy as np
import GMatElastic.Cartesian3d as GMat
import unittest

class Test(unittest.TestCase):

    def test_main(self):

        with h5py.File('Cartesian3d_random.hdf5', 'r') as data:

            mat = GMat.Array2d(data['/shape'][...])

            I = data['/model/I'][...]
            K = data['/model/K'][...]
            G = data['/model/G'][...]

            mat.setElastic(I, K, G)

            for i in range(20):

                GradU = data['/random/{0:d}/GradU'.format(i)][...]

                Eps = np.einsum('...ijkl,...lk->...ij', mat.I4s(), GradU)
                mat.setStrain(Eps)

                self.assertTrue(np.allclose(mat.Stress(), data['/random/{0:d}/Stress'.format(i)][...]))
                self.assertTrue(np.allclose(mat.Tangent(), data['/random/{0:d}/Tangent'.format(i)][...]))

if __name__ == '__main__':

    unittest.main()
