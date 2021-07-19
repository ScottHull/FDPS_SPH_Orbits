import os
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile

times = [0, 20, 1500, 3000]
number_processes = 100
path = "/scratch/shull4/gi_new"

fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
axs = [ax1, ax2, ax3, ax4]
fig.supxlabel('x')
fig.supylabel('y')

for index, time in enumerate(times):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    axs[index].scatter(
        [p.position[0] for p in particles if p.tag == 0],
        [p.position[1] for p in particles if p.tag == 0],
        marker="+",
        color='blue',
        label="Target Silicate"
    )
    axs[index].scatter(
        [p.position[0] for p in particles if p.tag == 2],
        [p.position[1] for p in particles if p.tag == 2],
        marker="+",
        color='green',
        label="Impactor Silicate"
    )
    axs[index].scatter(
        [p.position[0] for p in particles if p.tag == 1],
        [p.position[1] for p in particles if p.tag == 1],
        marker="+",
        color='red',
        label="Target Iron"
    )
    axs[index].scatter(
        [p.position[0] for p in particles if p.tag == 3],
        [p.position[1] for p in particles if p.tag == 3],
        marker="+",
        color='yellow',
        label="Impactor Iron"
    )
    axs[index].grid()
    axs[index].set_title("{} hrs".format(round(pm.time * 0.000277778, 4)))
    axs[index].set_xlim(-1e8, 1e8)
    axs[index].set_ylim(-1e8, 1e8)
    if index == 0:
        axs[index].legend(loc='lower right')

plt.savefig("multitime.png", format='png')

