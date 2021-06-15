import os
import shutil
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.geometry import get_impact_geometry, get_velocity_profile

start_time = 0
end_time = 100
interval = 5
number_processes = 100
path = "/scratch/shull4/gi"
output = "/scratch/shull4/impact_geometry"
velocity_output = "/scratch/shull4/velocity_impact_geometry"

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)
if os.path.exists(velocity_output):
    shutil.rmtree(velocity_output)
os.mkdir(velocity_output)

target_radius = 0
impactor_radius = 0
v_target_at_impact = 0  # desired value
v_impactor_at_impact = 0  # desired value

v_target_profile = []
v_impactor_profile = []
time_profile = []

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)
    
    if time == 0:
        G = 6.67 * 10 ** -11
        target = [p for p in particles if p.tag <= 1]
        impactor = [p for p in particles if p.tag > 1]
        target_radius = (max([p.distance for p in target]) - min([p.distance for p in target])) / 2.0
        impactor_radius = (max([p.distance for p in impactor]) - min([p.distance for p in impactor])) / 2.0

        tar_mass = sum([p.mass for p in target])
        imp_mass = sum([p.mass for p in impactor])
        tot_mass = tar_mass + imp_mass
        v_esc = sqrt((2 * G * tot_mass) / (target_radius + impactor_radius))

        v_target_at_impact = [
            (imp_mass / tot_mass) * v_esc,
            0,
            0
        ]
        v_impactor_at_impact = [
            (tar_mass / tot_mass) * v_esc,
            0,
            0
        ]

    target, impactor, target_com_x, target_com_y, target_com_z, impactor_com_x, impactor_com_y, impactor_com_z, \
    imp_angle, min_x_tar, max_x_tar, min_y_tar, max_y_tar, min_x_imp, max_x_imp, min_y_imp, max_y_imp = \
        get_impact_geometry(particles=particles)

    tar_velocity, imp_velocity, v_esc = get_velocity_profile(
        particles=particles,
        target_radius=target_radius,
        impactor_radius=impactor_radius
    )
    v_target_profile.append(tar_velocity)
    v_impactor_profile.append(imp_velocity)
    time_profile.append(time)

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
    ax.set_xlim(-2.1e7, 2.1e7)
    ax.set_ylim(-2.1e7, 2.1e7)
    ax.set_title("Iteration: {} // Impact Angle: {}".format(time, imp_angle))
    ax.grid()
    ax.legend(loc='upper left')
    plt.savefig(output + "/{}.png".format(time), format='png')



animate(
    start_time=0,
    end_time=end_time,
    interval=interval,
    path=output,
    fps=5,
    filename="impact_geometry.mp4",
)

fig = plt.figure(figsize=(16, 9))
ax_x = fig.add_subplot(311)
ax_y = fig.add_subplot(312)
ax_z = fig.add_subplot(313)
ax_x.plot(
    time_profile,
    [i[0] for i in v_target_profile],
    linewidth=2.0,
    color='blue',
    label="TARGET"
)
ax_x.plot(
    time_profile,
    [i[0] for i in v_impactor_profile],
    linewidth=2.0,
    color='red',
    label="IMPACTOR"
)
ax_x.axhline(
    v_target_at_impact[0],
    linewidth=2.0,
    linestyle="--",
    color='blue',
    label="EXPECTED V TARGET AT IMPACT"
)
ax_x.axhline(
    v_impactor_at_impact[0],
    linewidth=2.0,
    linestyle="--",
    color='red',
    label="EXPECTED V IMPACTOR AT IMPACT"
)

ax_y.plot(
    time_profile,
    [i[1] for i in v_target_profile],
    linewidth=2.0,
    color='blue',
    label="TARGET"
)
ax_y.plot(
    time_profile,
    [i[1] for i in v_impactor_profile],
    linewidth=2.0,
    color='red',
    label="IMPACTOR"
)
ax_y.axhline(
    v_target_at_impact[1],
    linewidth=2.0,
    linestyle="--",
    color='blue',
    label="EXPECTED V TARGET AT IMPACT"
)
ax_y.axhline(
    v_impactor_at_impact[1],
    linewidth=2.0,
    linestyle="--",
    color='red',
    label="EXPECTED V IMPACTOR AT IMPACT"
)

ax_z.plot(
    time_profile,
    [i[2] for i in v_target_profile],
    linewidth=2.0,
    color='blue',
    label="TARGET"
)
ax_z.plot(
    time_profile,
    [i[2] for i in v_impactor_profile],
    linewidth=2.0,
    color='red',
    label="IMPACTOR"
)
ax_z.axhline(
    v_target_at_impact[2],
    linewidth=2.0,
    linestyle="--",
    color='blue',
    label="EXPECTED V TARGET AT IMPACT"
)
ax_z.axhline(
    v_impactor_at_impact[2],
    linewidth=2.0,
    linestyle="--",
    color='red',
    label="EXPECTED V IMPACTOR AT IMPACT"
)

ax_z.set_xlabel("Iteration")
ax_x.set_ylabel("v_x")
ax_y.set_ylabel("v_y")
ax_z.set_ylabel("v_z")
ax_x.set_title("Impact Velocity Profile")
ax_x.grid()
ax_y.grid()
ax_z.grid()
ax_z.legend()

plt.savefig("impact_velocity.png", format='png')

