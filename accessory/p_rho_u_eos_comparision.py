import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("dark_background")


new_phase_path = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH_Orbits/src/phase_data/forst_STS.rho_u.txt"
old_phase_path = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH_Orbits/src/phase_data/duniteN.rho_u.txt"

new_phase_df = pd.read_fwf(new_phase_path, skiprows=2,
                        names=["density", "internal_energy", "temperature", "pressure", "soundspeed", "entropy"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=2,
                        names=["density", "internal_energy", "temperature", "pressure", "soundspeed", "entropy"])


fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    new_phase_df['internal_energy'],
    [p / rho for p, rho in zip(new_phase_df['pressure'], new_phase_df['density'])],
    s=2,
    c='blue',
    label="Stewart M-ANEOS"
)
ax.scatter(
    old_phase_df['internal_energy'],
    [p / rho for p, rho in zip(old_phase_df['pressure'], old_phase_df['density'])],
    s=2,
    c='red',
    label="GADGET M-ANEOS"
)
ax.set_xlabel("Internal Energy")
ax.set_ylabel(r"P/$\rho$")
ax.grid(alpha=0.4)
ax.legend()

plt.show()
