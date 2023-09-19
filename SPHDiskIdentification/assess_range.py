#!/usr/bin/env python3
import os
import shutil
import pandas as pd
from random import randint
import matplotlib.pyplot as plt

from src.combine import CombinedFile
from src.identify import ParticleMap
from src.animate import animate

"""
This file can be used to profile a range of timesteps for a FDPS SPH generated file.
"""

# define where the data is
path = '/home/theia/scotthull/Paper2_SPH/gi/500_mars/500_mars'
to_path = 'test_animation'
start_iteration = 60
end_iteration = 1500
increment = 20
number_of_processes = 600
file_headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                "potential energy", "entropy", "temperature"]

# define the planet parameters
mass_planet = 6.39e23
equatorial_radius = 3390e3

if os.path.exists(to_path):
    shutil.rmtree(to_path)
os.mkdir(to_path)

for iteration in range(start_iteration, end_iteration + 1, increment):
    print(f"Processing iteration {iteration}...")
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

    plt.savefig(f"{to_path}/{iteration}.png", dpi=300)

    print(f"Finished iteration {iteration}!\n")

animate(
    start_time=start_iteration,
    end_time=end_iteration,
    interval=increment,
    path=to_path,
    fps=15
)
