import GMatElastic as gmat
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.style.use(["goose", "goose-latex"])
except FileNotFoundError:
    pass


def ddot42(A4, B2):
    return np.einsum("ijkl,lk->ij", A4, B2)


def ddot22(A2, B2):
    return np.einsum("ij,ji", A2, B2)


I4d = gmat.Cartesian3d.I4d()

mat = gmat.Cartesian3d.Elastic(10.0, 1.0)

epseq = np.zeros(101)
sigeq = np.zeros(101)

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, len(epseq))):

    Eps = np.array(
        [
            [0.0, gamma, 0.0],
            [gamma, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    Sig = mat.Stress(Eps)

    Epsd = ddot42(I4d, Eps)
    Sigd = ddot42(I4d, Sig)

    epseq[igamma] = np.sqrt(2.0 / 3.0 * ddot22(Epsd, Epsd))
    sigeq[igamma] = np.sqrt(3.0 / 2.0 * ddot22(Sigd, Sigd))

# Plot

fig, ax = plt.subplots()

ax.plot(epseq, sigeq)

ax.set_xlabel(r"$\varepsilon_\mathrm{eq}$")
ax.set_ylabel(r"$\sigma_\mathrm{eq}$")

plt.savefig("stress-strain.pdf")
plt.show()
