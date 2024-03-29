import os
import csv
from math import atan, sin, pi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.animate import animate

# use dark background
# plt.style.use("dark_background")

base_path = "/Users/scotthull/Downloads/mars_data_robin/csv/"
headers = [  # tag: 1=water ; 2=ice ; 3=serpentine ; 4=dunite ; 5=iron
    'id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'temperature', 'size', 'tag'
]
square_scale = 100
fname_template = "marsimp{}.csv"

mars_radius = 3390  # km

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


def get_impact_angle(com_target, com_impactor):
    """
    Given the center of mass between the target and impactor as cartesian coordinates, calculate the angle between them.
    :param com_target:
    :param com_impactor:
    :return:
    """
    # calculate the angle between the target and impactor
    impact_angle = atan((com_impactor[1] - com_target[1]) / (com_impactor[0] - com_target[0]))
    # return impact parameter b = sin(impact_angle)
    return np.sin(impact_angle)

def get_target_impactor_separation(com_target: tuple, com_impactor: tuple):
    """
    Calculate the distance between the target and impactor COMs and the x and y components
    :param com_target:
    :param com_impactor:
    :return:
    """
    # calculate the distance between the target and impactor COMs
    x_sep = com_impactor[0] - com_target[0]
    y_sep = com_impactor[1] - com_target[1]
    total_distance = np.sqrt(x_sep ** 2 + y_sep ** 2)
    return x_sep, y_sep, total_distance

# the files in the base path have a number in the string
# we need to order the file names by this number
dir_files = [i for i in os.listdir(base_path) if i.endswith(".csv")]
dir_files.sort(key=lambda x: int(x.replace("marsimp", "").replace(".csv", "")))


def convert_fwf_to_csv(fwf_path: str):
    """
    Read in a FWF file using the CSV module and write out a CSV file.
    :param fwf_path:
    :return:
    """
    # read in the FWF file
    with open(fwf_path, 'r') as fwf_file:
        # read in the lines
        lines = fwf_file.readlines()
        # open the CSV file
        with open(fwf_path.replace(".csv", ".csv"), 'w') as csv_file:
            # create a CSV writer object
            csv_writer = csv.writer(csv_file)
            # loop through each line
            for line in lines:
                # write out the line as a CSV
                csv_writer.writerow(line.split())


# for index, i in enumerate(dir_files):
#     convert_fwf_to_csv(base_path + i)

target_particles = None
impactor_particles = None
target_velocity = []
impactor_velocity = []
impact_angle = []
iterations = 0
# loop through all the files in the directory
# for index, i in enumerate(dir_files):
for index, i in enumerate(np.arange(1, 5, 1)):
    i = fname_template.format(str(i).zfill(3))
    # read in as a dataframe
    df = pd.read_csv(f"{base_path}{i}", header=None)
    # set the column names
    df.columns = headers
    utime = 0.00025826343 ** -1
    df['mass'] = 1.0
    df['vx'] *= (10 ** 3) * (1 / utime)
    df['vy'] *= (10 ** 3) * (1 / utime)
    df['vz'] *= (10 ** 3) * (1 / utime)
    df['x'] *= (10 ** 3)
    df['y'] *= (10 ** 3)
    df['z'] *= (10 ** 3)
    df['velocity'] = np.sqrt(df['vx'] ** 2 + df['vy'] ** 2 + df['vz'] ** 2)
    df['radius'] = np.sqrt(df['x'] ** 2 + df['y'] ** 2 + df['z'] ** 2)
    if index == 0:
        target = df[df['x'] <= 3.5]
        impactor = df[df['x'] > 3.5]
        target_particles = target['id']
        impactor_particles = impactor['id']
    target = df[df['id'].isin(target_particles)]
    impactor = df[df['id'].isin(impactor_particles)]
    com_target = center_of_mass(
        target['x'].values, target['y'].values, target['z'].values, target['mass'].values
    )
    com_impactor = center_of_mass(
        impactor['x'].values, impactor['y'].values, impactor['z'].values, impactor['mass'].values
    )
    df['x'] -= com_target[0]
    df['y'] -= com_target[1]
    df['z'] -= com_target[2]
    target = df[df['id'].isin(target_particles)]
    impactor = df[df['id'].isin(impactor_particles)]
    com_target = center_of_mass(
        target['x'].values, target['y'].values, target['z'].values, target['mass'].values
    )
    com_impactor = center_of_mass(
        impactor['x'].values, impactor['y'].values, impactor['z'].values, impactor['mass'].values
    )
    num_target_particles = len(target)
    num_impactor_particles = len(impactor)
    total_particles = num_target_particles + num_impactor_particles
    gamma = num_impactor_particles / num_target_particles
    print(f"Gamma: {gamma}")
    impact_angle.append(get_impact_angle(center_of_mass(target['x'], target['y'], target['z'], target['mass']), center_of_mass(impactor['x'], impactor['y'], impactor['z'], impactor['mass'])))
    target_velocity.append(target['velocity'].mean())
    impactor_velocity.append(impactor['velocity'].mean())
    com_x, com_y, com_z = center_of_mass(
        df['x'].values, df['y'].values, df['z'].values, df['mass'].values
    )
    if index == 0:
        # calculate the inital distance between the target and impactor COMs
        x_sep, y_sep, total_distance = get_target_impactor_separation(com_target, com_impactor)
        print(f"Initial distance between target and impactor: {total_distance / mars_radius:.2f} M_Mars (x component: {x_sep / mars_radius:.2f} M_Mars, y component: {y_sep / mars_radius:.2f} M_Mars)")
        print("Initial target position ", com_target)
        print("Initial impactor position", com_impactor)
    # print(f"Center of mass: {com_x}, {com_y}, {com_z}")
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    for j in [target, impactor]:
        ax.scatter(
            j['x'] - com_x, j['y'] - com_y, marker=".", s=5
        )
    # scatter the target and impactor COMs
    ax.scatter(
        com_target[0] - com_x, com_target[1] - com_y, marker="*", s=100, label="COM Target"
    )
    ax.scatter(
        com_impactor[0] - com_x, com_impactor[1] - com_y, marker="*", s=100, label="COM Impactor"
    )
    # turn off the axes
    # ax.set_xticks([], minor=False)
    # ax.set_yticks([], minor=False)
    # ax.set_xlim(-square_scale, square_scale)
    # ax.set_ylim(-square_scale, square_scale)
    ax.set_title(i)

    legend = plt.legend(loc='upper right', fontsize=16, title="Tag")
    # make markers larger in the legend
    for handle in legend.legendHandles:
        handle.set_sizes([200.0])
    plt.tight_layout()
    # plt.savefig(f"{base_path}{index}.png", dpi=200)
    # plt.show()

    print("Target Velocity: {} m/s".format(target['velocity'].mean()) * 1000)
    print("Impactor Velocity: {} m/s".format(impactor['velocity'].mean()) * 1000)
    print("Impact parameter: {:.3f}".format(sin(get_impact_angle(com_target, com_impactor))))

    iterations += 1

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
for i, label in [(target_velocity, "Target"), (impactor_velocity, "Impactor")]:
    ax.plot(range(0, iterations), np.array(i), linewidth=1.0, label=label)
    # scatter the points on top and annotate the value
    ax.scatter(
        range(0, iterations), np.array(i)
    )
    for point in range(0, iterations):
        ax.annotate(
            f"{np.array(i)[point]:.2f}",
            (point, np.array(i)[point]),
            textcoords="offset points",
            xytext=(0, 10),
            ha='center'
        )
ax.set_xlabel("Iteration")
ax.set_ylabel("Velocity (km/s)")
ax.legend()
ax.grid()
plt.show()

# plot the impact angle as a function of time
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.plot(range(0, iterations), np.array(impact_angle), linewidth=1.0)
# scatter the points on top and annotate the value
ax.scatter(
    range(0, iterations), np.array(impact_angle)
)
for point in range(0, iterations):
    ax.annotate(
        f"{np.array(impact_angle)[point]:.2f}",
        (point, np.array(impact_angle)[point]),
        textcoords="offset points",
        xytext=(0, 10),
        ha='center'
    )
ax.set_xlabel("Iteration")
ax.set_ylabel("Impact Parameter (b)")
ax.set_ylim(0, 1)
ax.grid()
plt.show()

animate(
    0,
    len([i for i in os.listdir(base_path) if i.endswith(".csv")]),
    interval=1,
    path=base_path,
    fps=5,
    filename=f"{base_path}_robin_impact.mp4"
)
