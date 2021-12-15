import numpy as np
import matplotlib.pyplot as plt

plt.style.use("dark_background")

alpha = 10 ** -5
deltaT = 10 ** 4
rho_0 = 7500

def calculate_delta_rho(alpha, rho_0, deltaT):
    return alpha * rho_0 * deltaT

deltaT_range = np.arange(1000, 40000, 100)
alpha_range = np.arange(1e-6, 1e-4, 1e-7)

t_sens = [calculate_delta_rho(alpha=alpha, rho_0=rho_0, deltaT=t) for t in deltaT_range]
alpha_sens = [calculate_delta_rho(alpha=a, rho_0=rho_0, deltaT=deltaT) for a in alpha_range]

fig, axs = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.0, "wspace": 0.14})
ax1, ax2 = axs.flatten()
ax1.plot(deltaT_range, t_sens, linewidth=2.0)
ax2.plot(alpha_range, alpha_sens, linewidth=2.0)
ax1.set_xlabel("delta T"), ax1.set_ylabel("delta Rho")
ax2.set_xlabel("alpha"), ax2.set_ylabel("delta Rho")
for ax in axs.flatten():
    ax.grid(alpha=0.4)

plt.show()
