import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

start_time = 0
end_time = 3000
interval = 20
number_processes = 100
path = "/scratch/shull4/gi_impact5"
output = "/scratch/shull4/low_density_entropy_evolutoim"

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)

cf_end = CombineFile(num_processes=number_processes, time=end_time, output_path=path)
formatted_time_end = cf_end.sim_time
combined_file_end = cf_end.combine()
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm_end = ParticleMap(path=f, center=False, relative_velocity=False)
particles = pm_end.collect_particles()
pm_end.solve(particles=particles)
os.remove(f)

end = {}
for p in particles:
    end.update({p.particle_id: p.label})

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    disk = [p for p in particles if end[p.particle_id] == "DISK" and p.density < 1-0]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.density for p in disk],
        [p.entropy for p in disk],
        marker="+",
        color="blue",
    )

    ax.set_xlabel("Density")
    ax.set_ylabel("Entropy")
    ax.set_title("Time: {} sec (iteration: {})".format(formatted_time, time))
    ax.grid()
    ax.legend(loc="upper left")
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
