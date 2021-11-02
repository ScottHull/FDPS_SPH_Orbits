import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from random import randint

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

start_time = 0
end_time = 3000
interval = 20
number_processes = 100
path = "/scratch/shull4/gi_impact5"
output = "/scratch/shull4/low_density_entropy_evolution"

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)

cf_end = CombineFile(num_processes=number_processes, time=end_time, output_path=path)
formatted_time_end = cf_end.sim_time
combined_file_end = cf_end.combine()
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm_end = ParticleMap(path=f, center=False, relative_velocity=False)
particles = pm_end.collect_particles()
pm_end.solve(particles=particles, phase_path="src/phase_data/forstSTS__vapour_curve.txt")
os.remove(f)

end = {}
for p in particles:
    end.update({p.particle_id: p.label})

tracked_particles = {}
tracked_iterations = {}
disk = [p for p in particles if end[p.particle_id] == "DISK" and p.density < 10]

for i in range(0, 10):
    r = disk[randint(0, len(disk) - 1)]
    tracked_particles.update({r.particle_id: []})
    tracked_iterations.update({r.particle_id: []})

cmap = plt.get_cmap("viridis")
norm = plt.Normalize(0, 10000)

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    disk = [p for p in particles if end[p.particle_id] == "DISK"]


    tp = [p for p in particles if p.particle_id in tracked_particles.keys()]
    for i in tp:
        tracked_particles[i.particle_id].append(i.entropy)
        tracked_iterations[i.particle_id].append(time)

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.position[0] for p in disk],
        [p.position[1] for p in disk],
        marker="+",
        c=[p.entropy for p in disk]
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.ax.set_title("Entropy")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Time: {} sec (iteration: {})".format(formatted_time, time))
    ax.grid()
    ax.set_xlim(-1e8, 1e8)
    ax.set_ylim(-1e8, 1e8)

    plt.savefig(output + "/{}.png".format(time), format='png')

animate(
    start_time=0,
    end_time=end_time,
    interval=interval,
    path=output,
    fps=5,
    filename="low_density_entropy_evolution.mp4",
)

shutil.rmtree(output)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
for i in tracked_particles.keys():
    ax.plot(
        tracked_iterations[i],
        tracked_particles[i],
        linewidth=2.0,
        label=i
    )

ax.set_xlabel("Iteration")
ax.set_ylabel("Entropy")
ax.set_title("Disk Particle Entropy Time Evolution")
ax.grid()
ax.legend()

plt.savefig("disk_entropy_time_evolution.png", format='png')
