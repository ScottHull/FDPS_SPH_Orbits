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
path = "/scratch/shull4/gi"
t_output = "/scratch/shull4/t_output"
s_output = "/scratch/shull4/s_output"

for i in [t_output, s_output]:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.position[0] / 1000 for p in particles],
        [p.position[1] / 1000 for p in particles],
        marker="+",
        c=[p.entropy for p in particles],
        cmap='viridis'
    )
    plt.colorbar()
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_title("Entropy at Iteration: {}".format(time))
    ax.grid()
    plt.savefig(s_output + "/{}.png".format(time), format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.position[0] / 1000 for p in particles],
        [p.position[1] / 1000 for p in particles],
        marker="+",
        c=[p.temperature for p in particles],
        cmap='viridis'
    )
    plt.colorbar()
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_title("Temperature at Iteration: {}".format(time))
    ax.grid()
    plt.savefig(t_output + "/{}.png".format(time), format='png')
