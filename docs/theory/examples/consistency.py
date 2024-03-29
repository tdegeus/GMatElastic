import GMatElastic as gmat
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose", "goose-latex"])


def ddot42(A4, B2):
    return np.einsum("ijkl,lk->ij", A4, B2)


def ddot22(A2, B2):
    return np.einsum("ij,ji", A2, B2)


def norm(A2):
    return np.abs(np.einsum("ij,ji", A2, A2))


mat = gmat.Cartesian3d.Elastic(10.0, 1.0)

Eps_star = np.array(
    [
        [0.0, 0.1, 0.0],
        [0.1, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
)

Sig_star, C4 = mat.Tangent(Eps_star)

x = np.logspace(-16, 0, 100)
y = np.zeros(x.shape)

for i in range(len(x)):

    delta_Eps = np.random.random((3, 3)) * x[i]
    delta_Eps = 0.5 * (delta_Eps + delta_Eps.T)

    Sig = mat.Stress(Eps_star + delta_Eps)

    delta_Sig = Sig - Sig_star

    y[i] = norm(delta_Sig - ddot42(C4, delta_Eps)) / norm(delta_Sig)

# Plot

fig, ax = plt.subplots()

ax.plot(x, y, color="r", label=r"measurement")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$|| \delta \bm{\varepsilon} ||$")
ax.set_ylabel(r"$\eta$")

gplt.plot_powerlaw(
    -2,
    0,
    1,
    1,
    axis=ax,
    color="k",
    linewidth=1,
    label=r"rounding error: $|| \delta \bm{\varepsilon} ||^{-2}$",
)

ax.legend()

plt.savefig("consistency.pdf")
plt.show()
