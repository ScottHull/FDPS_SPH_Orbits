#!/usr/bin/env python3
import pandas as pd
from random import randint
import matplotlib.pyplot as plt

from src.combine import CombinedFile
from src.identify import ParticleMap

"""
This file can be used to profile a single timestep for a FDPS SPH generated file.
"""

# define where the data is
path = '/home/theia/scotthull/Paper2_SPH/gi/500_mars/500_mars'
number_of_processes = 600
iteration = 500
file_headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                "potential energy", "entropy", "temperature"]

# define the planet parameters
mass_planet = 6.39e23
equatorial_radius = 3390e3

# create the combined file
c = CombinedFile(
    path=path,
    iteration=iteration,
    number_of_processes=number_of_processes,
    to_fname=f"merged_{iteration}_{randint(1, int(1e5))}.dat"
)
combined_file = c.combine_to_memory()
# replace the headers
combined_file.columns = file_headers
time = c.sim_time
# create the particle map
particle_map = ParticleMap(particles=combined_file, mass_planet=mass_planet, equatorial_radius=equatorial_radius)
particles = particle_map.loop()

print("Particle map created!")
print(
    f"Time: {time} hours\n"
    f"Iteration: {iteration}\n"
    f"Planet mass: {particle_map.mass_planet} kg\n"
    f"Planet radius: {particle_map.equatorial_radius / 1000} km\n"
    f"Planet density: {particle_map.bulk_density} kg/m3\n"
    f"# escaping particles: {len(particles[particles['label'] == 'ESCAPE'])}\n"
    f"# disk particles: {len(particles[particles['label'] == 'DISK'])}\n"
    f"# planet particles: {len(particles[particles['label'] == 'PLANET'])}\n"
    f"# error particles: {len(particles[particles['label'] == None])}\n"
    f"Total particles: {len(particles)}\n"
    f"Target particles: {len(particles[particles['tag'] < 2])}\n"
    f"Impactor particles: {len(particles[particles['tag'] > 1])}\n"
    f"Target Mass: {sum(particles[particles['tag'] < 2]['mass'])} kg\n"
    f"Impactor Mass: {sum(particles[particles['tag'] > 1]['mass'])} kg\n"
)

fig = plt.figure(figsize=(10, 10))
# use the dark background
plt.style.use('dark_background')
ax = fig.add_subplot(111)
ax.set_title(f"{time} hrs. (Iteration {iteration})")
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
# scatter planet, disk, and escaping particles
for i in ['PLANET', 'DISK', 'ESCAPE']:
    ax.scatter(
        particles[particles['label'] == i]['x'] / 10 ** 7,
        particles[particles['label'] == i]['y'] / 10 ** 7,
        marker='.',
        s=2,
        label=i
    )
legend = ax.legend(loc='upper right')
# increase legend marker size
for i in legend.legendHandles:
    i.set_sizes([20])

plt.savefig(f"test_{iteration}.png", dpi=300)
