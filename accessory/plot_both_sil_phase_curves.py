import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("dark_background")

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    new_phase_df['entropy_sol_liq'],
    new_phase_df['temperature'],
    linewidth=2.0,
    label="New Liquid",
)
ax.plot(
    new_phase_df['entropy_vap'],
    new_phase_df['temperature'],
    linewidth=2.0,
    label="New Vapor"
)
ax.plot(
    old_phase_df['entropy_sol_liq'],
    old_phase_df['temperature'],
    linewidth=2.0,
    label="Old Liquid",
)
ax.plot(
    old_phase_df['entropy_vap'],
    old_phase_df['temperature'],
    linewidth=2.0,
    label="Old Vapor"
)
ax.set_xlim(-10, 15000)
ax.set_xlabel("Entropy")
ax.set_ylabel("Temperatuare")
ax.set_title("New vs. Old EoS Phase Curves")
ax.grid(alpha=0.4)
ax.legend()
plt.show()
