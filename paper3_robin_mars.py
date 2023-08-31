import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.animate import animate

# use dark background
plt.style.use("dark_background")

base_path = "/Users/scotthull/Downloads/data_robin/"
headers = [  # tag: 1=water ; 2=ice ; 3=serpentine ; 4=dunite ; 5=iron
    'id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'temperature', 'size', 'tag'
]
square_scale = 1e6

def center_of_mass(x: np.array, y: np.array, z: np.array, mass: np.array):
    """
    Calculate the center of mass of a system of particles
    :param x:
    :param y:
    :param z:
    :param mass:
    :return:
    """
    # calculate the center of mass
    x_com = np.sum(x * mass) / np.sum(mass)
    y_com = np.sum(y * mass) / np.sum(mass)
    z_com = np.sum(z * mass) / np.sum(mass)
    return x_com, y_com, z_com

# loop through all the files in the directory
for index, i in enumerate([i for i in os.listdir(base_path) if not i.endswith(".png")]):
    # read in as a dataframe
    df = pd.read_fwf(f"{base_path}{i}", header=None)
    # set the column names
    df.columns = headers
    df['mass'] = 1.0
    com_x, com_y, com_z = center_of_mass(
        df['x'].values, df['y'].values, df['z'].values, df['mass'].values
    )
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    for tag in df['tag'].unique():
        ax.scatter(
            df[df['tag'] == tag]['x'] - com_x, df[df['tag'] == tag]['y'] - com_y, marker=".", s=5, label=tag
        )
    # turn off the axes
    ax.set_xticks([], minor=False)
    ax.set_yticks([], minor=False)
    # ax.set_xlim(-square_scale, square_scale)
    # ax.set_ylim(-square_scale, square_scale)
    ax.set_title(i)

    legend = plt.legend(loc='upper right', fontsize=16, title="Tag")
    # make markers larger in the legend
    for handle in legend.legendHandles:
        handle.set_sizes([200.0])
    plt.tight_layout()
    plt.savefig(f"{base_path}{index}.png", dpi=200)

animate(
    0,
    len([i for i in os.listdir(base_path) if not i.endswith(".png")]),
    interval=1,
    path=base_path,
    fps=5,
    filename=f"{base_path}robin_impact.mp4"
)
