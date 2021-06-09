import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.geometry import get_impact_geometry

start_time = 0
end_time = 100
interval = 5
number_processes = 100
path = "/scratch/shull4/gi"
output = "/scratch/shull4/impact_geometry"

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    target, impactor, target_com_x, target_com_y, target_com_z, impactor_com_x, impactor_com_y, impactor_com_z, \
    imp_angle, min_x_tar, max_x_tar, min_y_tar, max_y_tar, min_x_imp, max_x_imp, min_y_imp, max_y_imp = \
        get_impact_geometry(particles=particles)

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.position[0] for p in target],
        [p.position[1] for p in target],
        marker="+",
        color='blue',
        label="TARGET"
    )
    ax.scatter(
        [p.position[0] for p in impactor],
        [p.position[1] for p in impactor],
        marker="+",
        color='red',
        label="IMPACTOR"
    )
    ax.plot(
        [target_com_x, impactor_com_x],
        [target_com_y, target_com_y],
        linewidth=2.0,
        color='black'
    )
    ax.plot(
        [impactor_com_x, impactor_com_x],
        [target_com_y, impactor_com_y],
        linewidth=2.0,
        color='black'
    )
    ax.plot(
        [target_com_x, impactor_com_x],
        [target_com_y, impactor_com_y],
        linewidth=2.0,
        color='black',
        label="IMP ANGLE: {}".format(imp_angle)
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y"),
    ax.set_title("Iteration: {} // Impact Angle: {}".format(time, imp_angle))
    ax.grid()
    ax.legend(loc='upper left')
    plt.savefig(output + "{}.png".format(time), format='png')

animate(
    start_time=0,
    end_time=end_time,
    interval=interval,
    path=output,
    fps=5,
    filename="impact_geometry.py.mp4",
)
