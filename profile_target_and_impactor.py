from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt

target_path = "/home/theia/scotthull/1M_high_rho_cutoff/tar.dat"
impactor_path = "/home/theia/scotthull/1M_high_rho_cutoff/imp.dat"
target_df = pd.read_csv(target_path, skiprows=2, header=None, delimiter="\t")
impactor_df = pd.read_csv(impactor_path, skiprows=2, header=None, delimiter="\t")
target_radius = [sqrt(i ** 2 + j ** 2 + k ** 2) / (6371 * 1000) for i, j, k in
                 zip(target_df[3], target_df[4], target_df[5])]
impactor_radius = [sqrt(i ** 2 + j ** 2 + k ** 2) / (6371 * 1000) for i, j, k in
                   zip(impactor_df[3], impactor_df[4], impactor_df[5])]

fig, axs = plt.subplots(2, 3, figsize=(16, 20), sharex='all', gridspec_kw={"hspace": 0.10, "wspace": 0.14})

ax1, ax2, ax3, ax4, ax5, ax6 = axs.flatten()
ax1.scatter(
    target_radius, target_df['entropy'], s=1
)
ax2.scatter(
    target_radius, target_df['temperature'], s=1
)
ax3.scatter(
    target_radius, target_df['density'], s=1
)
ax4.scatter(
    impactor_radius, impactor_df['entropy'], s=1
)
ax5.scatter(
    impactor_radius, impactor_df['temperature'], s=1
)
ax6.scatter(
    impactor_radius, impactor_df['density'], s=1
)
ax1.set_ylabel("Entropy (Target)"), ax2.set_ylabel("Temperature (Target)"), ax3.set_ylabel("Density (Target)")
ax4.set_ylabel("Entropy (Impactor)"), ax5.set_ylabel("Temperature (Impactor)"), ax6.set_ylabel("Density (Impactor)")
ax4.set_xlabel(r'Radius $R_{\bigoplus}$'), ax5.set_xlabel(r'Radius $R_{\bigoplus}$'), ax6.set_xlabel(
    r'Radius $R_{\bigoplus}$')
for ax in axs.flatten():
    ax.grid(alpha=0.4)

plt.savefig("profile_target_and_impactor.png", format='png')
