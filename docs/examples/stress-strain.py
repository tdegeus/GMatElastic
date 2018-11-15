import GMatLinearElastic as gmat
import matplotlib.pyplot as plt
import numpy             as np

plt.style.use(['goose', 'goose-latex'])

ddot42 = lambda A4,B2: np.einsum('ijkl,lk->ij',A4,B2)
ddot22 = lambda A2,B2: np.einsum('ij,ji',A2,B2)

I4d = gmat.Cartesian3d.I4d()

mat = gmat.Cartesian3d.Elastic(10., 1.)

epseq = np.zeros(101)
sigeq = np.zeros(101)

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, len(epseq))):

  Eps = np.array([
    [  0.0, gamma,   0.0],
    [gamma,   0.0,   0.0],
    [  0.0,   0.0,   0.0],
  ])

  Sig = mat.Sig(Eps)

  Epsd = ddot42(I4d,Eps)
  Sigd = ddot42(I4d,Sig)

  epseq[igamma] = np.sqrt(2./3.*ddot22(Epsd,Epsd))
  sigeq[igamma] = np.sqrt(3./2.*ddot22(Sigd,Sigd))

fig, ax = plt.subplots()

ax.plot(epseq, sigeq)

ax.set_xlabel(r'$\varepsilon_\mathrm{eq}$')
ax.set_ylabel(r'$\sigma_\mathrm{eq}$')

plt.savefig('stress-strain.pdf')
plt.show()
