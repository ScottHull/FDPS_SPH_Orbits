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
path = "/scratch/shull4/gi_new"
output = "/scratch/shull4/map_final_to_all"

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

    planet = [p for p in particles if end[p.particle_id] == "PLANET"]
    disk = [p for p in particles if end[p.particle_id] == "DISK"]
    escape = [p for p in particles if end[p.particle_id] == "ESCAPE"]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.position[0] for p in planet],
        [p.position[1] for p in planet],
        marker="+",
        color="blue",
        label="PLANET"
    )
    ax.scatter(
        [p.position[0] for p in disk],
        [p.position[1] for p in disk],
        marker="+",
        color="green",
        label="DISK"
    )
    ax.scatter(
        [p.position[0] for p in escape],
        [p.position[1] for p in escape],
        marker="+",
        color="red",
        label="ESCAPE"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
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
    filename="map_final_to_all.mp4",
)
