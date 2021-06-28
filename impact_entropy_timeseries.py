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
end_time = 30
interval = 1
number_processes = 100
path = "/scratch/shull4/gi_impact5"

times = []
averages = []

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    target_silicate = [p.entropy for p in particles if p.tag == 0]
    target_iron = [p.entropy for p in particles if p.tag == 1]
    impactor_silicate = [p.entropy for p in particles if p.tag == 2]
    impactor_iron = [p.entropy for p in particles if p.tag == 3]

    times.append(time)
    averages.append(
        [
            sum(target_silicate) / len(target_silicate),
            sum(target_iron) / len(target_iron),
            sum(impactor_silicate) / len(impactor_silicate),
            sum(impactor_iron) / len(impactor_iron),
        ]
    )

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    times,
    [i[0] for i in averages],
    linewidth=2.0,
    label="Target Silicate"
)
ax.plot(
    times,
    [i[1] for i in averages],
    linewidth=2.0,
    label="Target Iron"
)
ax.plot(
    times,
    [i[2] for i in averages],
    linewidth=2.0,
    label="Impactor Silicate"
)
ax.plot(
    times,
    [i[3] for i in averages],
    linewidth=2.0,
    label="Impactor Iron"
)
ax.set_xlabel("Iteration")
ax.set_ylabel("Entropy")
ax.grid()
ax.legend()

plt.savefig("entropy_impact_evolution.png", format='png')
