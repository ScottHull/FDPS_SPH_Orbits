import numpy as np
from math import pi, log10
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')


target_mass = 5.29436E+24
impactor_mass = 7.9111E+23
densities = np.arange(1, 2000 + 1, 1)
particle_counts = [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8]

def get_particle_mass(body_mass, num_particles):
    return body_mass / num_particles

def get_smth_length(m, rho, xi=1.2, n=3):
    """
    h = xi * (m_j / rho_j)^(1/n)
    where xi is a coefficient, m_j is particle mass, and n is the resolution of the simulation.
    Increasing resolution decreases particle mass and therefore decreases smoothing length.
    Decreasing cutoff density increases smoothing length.
    """
    return xi * ((m / rho) ** (1 / n))

def wendland_c6(r, h, support_radius=2.5, n=3):
    h *= support_radius  # add support radius
    sigma = 1 / (h ** n)
    s = abs(r) / h
    W = (1365 * sigma) / (64 * pi)
    if s < 1:
        W *= ((1 - s) ** 8) * (1 + (8 * s) + (25 * (s ** 2)) + (32 * (s ** 3)))
    else:
        W *= 0
    return W



fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.set_xlabel("Density (kg/m$^3$)")
ax.set_ylabel("$h$ (km)")
# ax.set_title("Smoothing Length for Particle Counts (N)")
print(get_smth_length(get_particle_mass(target_mass + impactor_mass, 1.1e6), 5) / 1000)
print(get_smth_length(get_particle_mass(target_mass + impactor_mass, 1.1e6), 2000) / 1000)
for n in particle_counts:
    m = get_particle_mass(target_mass + impactor_mass, n)
    smth = [get_smth_length(m, rho) / 1000 for rho in densities]
    ax.plot(
        densities, smth, linewidth=2.0, label="N = %.1E" % n
    )
    ax.grid(alpha=0.4)
ax.legend()

# plt.show()
plt.savefig("smth_length.png", format='png', dpi=200)

# fig = plt.figure(figsize=(16, 9))
# ax = fig.add_subplot(111)
# ax.set_xlabel("Density (kg/$m^3$)")
# ax.set_ylabel("log(W(r=1, h))")
# ax.set_title("Wendland C6 Kernel for Particle Counts (N) and r=1")
# for n in particle_counts:
#     m = get_particle_mass(target_mass + impactor_mass, n)
#     smth = [get_smth_length(m, rho) for rho in densities]
#     wendland = [log10(wendland_c6(1, s)) for s in smth]
#     ax.plot(
#         densities, wendland, linewidth=2.0, label="N = %.1E" % n
#     )
#     ax.grid(alpha=0.4)
# ax.legend()

# plt.show()
# plt.savefig("wendlandc6.png", format='png', dpi=200)

# what happens when r_i - r_j -> \infty?
# lets assume that h = 100km constant
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.set_xlabel(r"$r_{i} - r_{j}$ (km)")
ax.set_ylabel("W($r_{i} - r_{j}$, $h$)")
# ax.set_title(r"Wendland C6 Kernel for Particle Counts (N)")
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for cd in [5, 2000]:
    if cd == 5:
        linestyle = "-"
    else:
        linestyle = "--"
    for n in particle_counts:
        particle_distances = [10**0, 10**1, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10]
        m = get_particle_mass(target_mass + impactor_mass, n)
        smth = get_smth_length(m, cd)
        wendland = [wendland_c6(r_i_minus_r_j, smth) for r_i_minus_r_j in particle_distances]
        ax.plot(
            [i / 1000 for i in particle_distances], wendland, linewidth=2.0, color=colors[particle_counts.index(n)], linestyle=linestyle
        )
        ax.grid(alpha=0.4)
for n in particle_counts:
    ax.scatter([], [], marker='s', s=80, color=colors[particle_counts.index(n)], label="N = %.1E" % n)
ax.plot([], [], linewidth=2.0, linestyle="-", color='black', label=r"${\rho}_{c}$ = 5 kg/m$^3$")
ax.plot([], [], linewidth=2.0, linestyle="--", color='black', label=r"${\rho}_{c}$ = 2000 kg/m$^3$")
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

plt.savefig("kernel.png", format='png', dpi=200)
