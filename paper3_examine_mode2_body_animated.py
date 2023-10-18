#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

from src.combine import CombineFile
from src.animate import animate

path = "/home/theia/scotthull/Paper3_SPH/tar-imp/target_mars_n_sph"
to_folder = "paper3_examine_mode2_body_animated"
start_iteration = 0
stop_iteration = 140
increment = 10
num_processes = 600
square_scale = 6 * 1e6

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]

if not os.path.exists(to_folder):
    os.mkdir(to_folder)

for iteration in np.arange(start_iteration, stop_iteration + increment, increment):
    cf = CombineFile(
        num_processes=num_processes, time=iteration, output_path=path, to_fname=""
    )
    df = cf.combine_df()
    df.columns = headers
    formatted_time = round(cf.sim_time * 0.000277778, 2)  # seconds -> hours
    # df['radius'] = np.sqrt(df['x'] ** 2 + df['y'] ** 2 + df['z'] ** 2)

    # first, make a 3D plot of the body in using xyz coordinates
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(df['x'], df['y'], df['z'], marker=".", s=5, alpha=1)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_zlim(-square_scale, square_scale)
    ax.set_title("Mars (t={} hrs.)".format(formatted_time))
    plt.tight_layout()
    plt.savefig(f"{to_folder}/{iteration}", dpi=150)

animate(
    start_time=start_iteration,
    end_time=stop_iteration,
    interval=increment,
    path=to_folder,
    fps=10,
    filename="paper3_examine_mode2_body_animated.mp4",
)
