import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn as sns

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.cross_section import cross_section_xy, sort_particles_by_closest


def format_seaborn_dataframe(old_eos_particles, new_eos_particles):
    particles = old_eos_particles + new_eos_particles
    new_or_old = ["old" for i in old_eos_particles] + ["new" for i in new_eos_particles]
    return pd.DataFrame({
        "particle_id": [p.particle_id for p in particles],
        "tag": [p.tag for p in particles],
        "new_or_old_eos": new_or_old,
        "label": [p.label for p in particles],
        "distance": [p.distance for p in particles],
        "entropy": [p.entropy for p in particles],
        "temperature": [p.temperature for p in particles],
        "pressure": [p.pressure for p in particles],
        "internal_energy": [p.internal_energy for p in particles],
        "density": [p.density for p in particles],
    })



new_gi_f= "/Users/scotthull/Desktop/gi_new_eos_3000.dat"
old_gi_f = "/Users/scotthull/Desktop/gi_old_eos_3000.dat"
pm_new_gi = ParticleMap(path=new_gi_f, center=True, relative_velocity=False)
pm_old_gi = ParticleMap(path=old_gi_f, center=True, relative_velocity=False)
particles_new_gi = pm_new_gi.collect_particles()
particles_old_gi = pm_old_gi.collect_particles()
# pm_end.solve(particles=particles)
# os.remove(f)

# particles = cross_section_xy(particles=particles, min_z=-100 * 10 ** 3, max_z=100 * 10 ** 3)
# particles = sort_particles_by_closest(particles=particles)

iron_particles_new_gi = [p for p in particles_new_gi if p.tag % 2 != 0 and p.distance < 1e7]
iron_particles_old_gi = [p for p in particles_old_gi if p.tag % 2 != 0 and p.distance < 1e7]
silicate_particles_new_gi = [p for p in particles_new_gi if p.tag % 2 == 0 and p.distance < 1e7]
silicate_particles_old_gi = [p for p in particles_old_gi if p.tag % 2 == 0 and p.distance < 1e7]


sns.set()

# fig, ax = plt.subplots(figsize=(16, 9))
# a0 = ax.scatter(
#     [p.position[0] for p in silicate_particles],
#     [p.position[1] for p in silicate_particles],
#     c=[p.entropy for p in silicate_particles],
#     cmap='Reds',
#     marker='o',
#     label="Silicate"
# )
# a1 = ax.scatter(
#     [p.position[0] for p in iron_particles],
#     [p.position[1] for p in iron_particles],
#     c=[p.entropy for p in iron_particles],
#     cmap='Blues',
#     marker='o',
#     label="Iron"
# )
# cbar1 = fig.colorbar(a0, ax=ax)
# cbar1.ax.set_ylabel('Entropy (Silicate)', rotation=270)
# cbar2 = fig.colorbar(a1, ax=ax)
# cbar2.ax.set_ylabel('Entropy (Iron)', rotation=270)
# ax.legend(loc='upper left')
# ax.set_xlim(-1e7, 1e7)
# ax.set_ylim(-1e7, 1e7)


df_iron = format_seaborn_dataframe(new_eos_particles=iron_particles_new_gi, old_eos_particles=iron_particles_old_gi)
df_silicate = format_seaborn_dataframe(new_eos_particles=silicate_particles_new_gi, old_eos_particles=silicate_particles_old_gi)


p = sns.jointplot(x='temperature', y='entropy', data=df_iron,
              hue="new_or_old_eos")
p.set_axis_labels('Temperature', 'Entropy', fontsize=14)
p.fig.suptitle("Iron")
p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95)  # Reduce plot to make room

p = sns.jointplot(x='temperature', y='entropy', data=df_silicate,
              hue="new_or_old_eos")
p.set_axis_labels('Temperature', 'Entropy', fontsize=14)
p.fig.suptitle("Silicate")
p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95)  # Reduce plot to make room

# ax = sns.kdeplot(x='distance', y='entropy', data=df, hue="tag", fill=True, cut=0, clip=(0, 1e7))



plt.show()
