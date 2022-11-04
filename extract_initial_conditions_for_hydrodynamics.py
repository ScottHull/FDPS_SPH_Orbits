#!/usr/bin/env python3
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt

runs = [
    "500_b073_new"
]

min_iteration = 0
max_iteration = 500

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
to_dir = "ic_for_hydrodynamics"
fig_dir = to_dir + "/figs"
if not os.path.exists(to_dir):
    os.mkdir(to_dir)

def get_time(f, local=True):
    formatted_time = None
    if local:  # if not reading from remote server
        with open(f, 'r') as infile:
            reader = csv.reader(infile, delimiter="\t")
            formatted_time = float(next(reader)[0])
        infile.close()
    else:
        formatted_time = float(next(f))
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours

for run in runs:
    if os.path.exists(to_dir + "/" + f"{run}.csv"):
        os.remove(to_dir + "/" + f"{run}")
    outfile = open(to_dir + "/" + run, 'w')
    writer = csv.writer(outfile, delimiter=",")
    writer.writerow(["run", "iteration", "time", "mean_velocity", "mean_density", "mean_entropy", "mean_internal_energy", "mean_pressure",
                     "mean_temperature"])

    path = base_path + f"/{run}formatted_{run}"
    # loop through all files in path, where the iteration is the file name minus the extension
    for file in os.listdir(path):
        f = path + "/" + file
        iteration = file.split(".")[0]
        if int(iteration) >= min_iteration and int(iteration) <= max_iteration:
            time = get_time(f)
            df = pd.read_csv(f, skiprows=2, sep="\t")
            # get the initial conditions for the hydrodynamics
            disk = df[df['label'] == "DISK"]
            mean_vel = disk['velocity'].mean()
            mean_density = disk['density'].mean()
            mean_entropy = disk['entropy'].mean()
            mean_internal_energy = disk['internal_energy'].mean()
            mean_pressure = disk['pressure'].mean()
            mean_temperature = disk['temperature'].mean()
            writer.writerow([run, time, mean_vel, mean_density, mean_entropy, mean_internal_energy, mean_pressure, mean_temperature])
    outfile.close()

# identify the max entropy of each run and write to a file, which you can use to identify the time of impact
outfile = open(to_dir + "/time_of_impact.csv", 'w')
writer = csv.writer(outfile, delimiter=",")
writer.writerow(["run", "iteration", "time", "mean_velocity", "mean_density", "mean_entropy", "mean_internal_energy", "mean_pressure",
                 "mean_temperature"])
for run in runs:
    df = pd.read_csv(to_dir + "/" + run, sep=",")
    max_entropy = df['mean_entropy'].max()
    time_of_impact = df[df['mean_entropy'] == max_entropy]['time'].values[0]
    # get the time of impact row as a list
    row = df[df['time'] == time_of_impact].values[0].tolist()
    # write row to file
    writer.writerow(row)

    # plot the time of impact
    # get the iteration of the time of impact
    iteration = df[df['time'] == time_of_impact]['iteration'].values[0]
    # get the file name of the time of impact
    file = base_path + f"/{run}formatted_{run}/" + f"{iteration}.csv"
    # read the file
    df = pd.read_csv(file, skiprows=2, sep="\t")
    planet, disk , escape = df[df['label'] == "PLANET"], df[df['label'] == "DISK"], df[df['label'] == "ESCAPE"]
    # use dark background
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(planet['x'], planet['y'], s=1,  label="planet")
    ax.scatter(disk['x'], disk['y'], s=1, label="disk")
    ax.scatter(escape['x'], escape['y'], s=1, label="escape")
    # title should include run, time, and iteration
    ax.set_title(f"{run} - time: {time_of_impact} hours, iteration: {iteration}")
    # save the figure as the run name
    plt.savefig(fig_dir + "/" + run + ".png", dpi=200)
outfile.close()
