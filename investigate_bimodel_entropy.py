#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.new_and_old_eos import seconds_to_hours

start_time = 0
end_time = 3000
interval = 20
number_processes = 200
path = "/home/theia/scotthull/gi_new_eos"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/track_high_entropy_particles"
density_output = "/home/theia/scotthull/FDPS_SPH_Orbits/track_high_entropy_particles_density"
density_output2 = "/home/theia/scotthull/FDPS_SPH_Orbits/track_high_entropy_particles_density2"

for o in [output, density_output]:
    if os.path.exists(o):
        shutil.rmtree(o)
    os.mkdir(o)

# plt.style.use("dark_background")
cf_end = CombineFile(num_processes=number_processes, time=end_time, output_path=path)
combined_file_end = cf_end.combine()
formatted_time_end = cf_end.sim_time
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm_end = ParticleMap(path=f, center=False, relative_velocity=False)
particles = pm_end.collect_particles()
pm_end.solve(particles=particles, phase_path="src/phase_data/forstSTS__vapour_curve.txt")
os.remove(f)

end = {}
high_entropy = {}
for p in particles:
    end.update({p.particle_id: p.label})
    if p.entropy > 8000 and p.label == "DISK":
        high_entropy.update({p.particle_id: p.entropy})
high_entropy_ids = list(high_entropy.keys())

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    planet = [p for p in particles if end[p.particle_id] == "PLANET"]
    disk = [p for p in particles if end[p.particle_id] == "DISK"]
    escape = [p for p in particles if end[p.particle_id] == "ESCAPE"]
    high_entropy_time = [p for p in particles if p.particle_id in high_entropy_ids]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    # fig.patch.set_facecolor('xkcd:black')
    ax.scatter(
        [p.position[0] for p in planet],
        [p.position[1] for p in planet],
        marker="o",
        color="blue",
        label="PLANET",
        s=0.2,
        alpha=0.4
    )
    ax.scatter(
        [p.position[0] for p in disk],
        [p.position[1] for p in disk],
        marker="o",
        color="green",
        label="DISK",
        s=0.2,
        alpha=0.4
    )
    ax.scatter(
        [p.position[0] for p in escape],
        [p.position[1] for p in escape],
        marker="o",
        color="red",
        label="ESCAPE",
        s=0.2,
        alpha=0.4
    )
    ax.scatter(
        [p.position[0] for p in high_entropy_time],
        [p.position[1] for p in high_entropy_time],
        marker="o",
        color="pink",
        label="S > 8000",
        s=2,
        alpha=1
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Time: {} sec (iteration: {})".format(round(seconds_to_hours(formatted_time), 2), time))
    ax.grid()
    ax.legend(loc="upper left")
    ax.set_xlim(-1e8, 1e8)
    ax.set_ylim(-1e8, 1e8)

    plt.savefig(output + "/{}.png".format(time), format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.density for p in disk],
        [p.entropy for p in disk],
        alpha=0.4,
        color='red',
        label="all disk particles"
    )
    ax.scatter(
        [p.density for p in high_entropy_time],
        [p.entropy for p in high_entropy_time],
        alpha=1,
        color='black',
        label="entropy > 8000"
    )
    ax.set_xlabel("Density")
    ax.set_ylabel("Entropy")
    ax.set_title("High Entropy Disk Particles")
    ax.grid()
    ax.legend()
    plt.savefig(density_output + "/{}.png".format(time), format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.density for p in disk if p.tag == 0],
        [p.entropy for p in disk if p.tag == 0],
        alpha=0.9,
        label="target silicate"
    )
    ax.scatter(
        [p.density for p in disk if p.tag == 1],
        [p.entropy for p in disk if p.tag == 1],
        alpha=0.9,
        label="target iron"
    )
    ax.scatter(
        [p.density for p in disk if p.tag == 2],
        [p.entropy for p in disk if p.tag == 2],
        alpha=0.9,
        label="impactor silicate"
    )
    ax.scatter(
        [p.density for p in disk if p.tag == 3],
        [p.entropy for p in disk if p.tag == 3],
        alpha=0.9,
        label="impactor iron"
    )
    ax.set_xlabel("Density")
    ax.set_ylabel("Entropy")
    ax.set_title("Disk Particles")
    ax.grid()
    ax.legend()
    plt.savefig(density_output2 + "/{}.png".format(time), format='png')

animate(
    start_time=start_time,
    end_time=end_time,
    interval=interval,
    path=output,
    fps=15,
    filename="high_entropy_evolution.mp4",
)

animate(
    start_time=start_time,
    end_time=end_time,
    interval=interval,
    path=density_output,
    fps=15,
    filename="density_high_entropy.mp4",
)

animate(
    start_time=start_time,
    end_time=end_time,
    interval=interval,
    path=density_output,
    fps=15,
    filename="density2_high_entropy.mp4",
)
