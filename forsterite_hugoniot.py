from src.hugoniot import Hugoniot_ANEOS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

forsterite_aneos_path = "src/phase_data/forst_STS.rho_u.txt"
forsterite_aneos_hugoniot_path = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH_Orbits/src/phase_data/forstSTS__hugoniot.txt"

h = Hugoniot_ANEOS()
h.read_ANEOS(forsterite_aneos_path)
rho1, P1, C1, U1, S1 = h.initial_conditions(forsterite_aneos_hugoniot_path)
T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s = h.rankine_hugoniot_equations(rho1, P1, U1)
P_s = np.array(P_s) / 10 ** 9
Us_s = np.array(Us_s) / 10 ** 3
Up_s = np.array(Up_s) / 10 ** 3
U_s = np.array(U_s) / 10 ** 6

fig, axs = plt.subplots(3, 2, figsize=(10, 15))
axs = axs.flatten()

axs[0].plot(
    Up_s, Us_s, linewidth=2.0, label="Calculated"
)
axs[0].plot(
    h.Up_h / 1000, h.Us_h / 1000, linewidth=2.0, linestyle="--", label="ANEOS"
)
axs[0].set_xlabel("Up (km/s)")
axs[0].set_ylabel("Us (km/s)")
axs[1].plot(
    P_s, T_s, linewidth=2.0, label="Calculated"
)
axs[1].plot(
    h.P_h / 10 ** 9, h.T_h, linewidth=2.0, linestyle="--", label="ANEOS"
)
axs[1].set_xlabel("P (GPa)")
axs[1].set_ylabel("T (K)")
axs[1].set_xlim(10**0, 10**4)
axs[2].plot(
    P_s, Rho_s, linewidth=2.0, label="Calculated"
)
axs[2].plot(
    h.P_h / 10 ** 9, h.rho_h, linewidth=2.0, linestyle="--", label="ANEOS"
)
axs[2].set_xlabel("P (GPa)")
axs[2].set_ylabel("Density (kg/m^3)")
axs[2].set_xlim(10**0, 10**4)
axs[3].plot(
    P_s, U_s, linewidth=2.0, label="Calculated"
)
axs[3].plot(
    h.P_h / 10 ** 9, h.U_h / 10 ** 6, linewidth=2.0, linestyle="--", label="ANEOS"
)
axs[3].set_xlabel("P (GPa)")
axs[3].set_ylabel("Internal Energy (J/K/kg)")
axs[3].set_xlim(10**0, 10**4)
axs[4].plot(
    Rho_s, T_s, linewidth=2.0, label="Calculated"
)
axs[4].plot(
    h.rho_h, h.T_h, linewidth=2.0, linestyle="--", label="ANEOS"
)
axs[4].set_xlabel("Density (kg/m^3)")
axs[4].set_ylabel("T (K)")

for ax in axs:
    ax.grid(alpha=0.4)
for ax in axs[1:-2]:
    ax.set_xscale("log")
for ax in axs[1:3]:
    ax.set_yscale("log")
axs[0].legend()

plt.show()



