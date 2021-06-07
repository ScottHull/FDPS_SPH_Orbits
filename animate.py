import os
import shutil
import numpy as np
from random import randint
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

start_time = 0
end_time = 3000
interval = 20
number_processes = 100
path = "/scratch/shull4/gi_v_esc"
output = "/scratch/shull4/animate"

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=True).collect_particles()
    os.remove(f)

    target_silicate = [p for p in pm if p.tag == 0]
    target_iron = [p for p in pm if p.tag == 1]
    impactor_silicate = [p for p in pm if p.tag == 2]
    impactor_iron = [p for p in pm if p.tag == 3]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.position[0] for p in target_silicate],
        [p.position[1] for p in target_silicate],
        marker="+",
        color='red',
        label="TARGET SILICATE"
    )
    ax.scatter(
        [p.position[0] for p in target_iron],
        [p.position[1] for p in target_iron],
        marker="+",
        color='blue',
        label="TARGET IRON"
    )
    ax.scatter(
        [p.position[0] for p in impactor_silicate],
        [p.position[1] for p in impactor_silicate],
        marker="+",
        color='purple',
        label="IMPACTOR SILICATE"
    )
    ax.scatter(
        [p.position[0] for p in impactor_iron],
        [p.position[1] for p in impactor_iron],
        marker="+",
        color='green',
        label="IMPACTOR IRON"
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
    filename="test_animate.mp4",
)
