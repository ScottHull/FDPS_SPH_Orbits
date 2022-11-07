#!/usr/bin/env python3
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt

runs = ['500_b073_new']

min_iteration = 0
max_iteration = 100

base_path = "/home/theia/scotthull/Paper1_SPH/gi"
to_dir = "ic_for_hydrodynamics"
fig_dir = to_dir + "/figs"

for d in [to_dir, fig_dir]:
    if not os.path.exists(d):
        os.mkdir(d)


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


def get_ic():
    outfile = open(to_dir + "/disk_history.csv", 'w')
    writer = csv.writer(outfile, delimiter=",")
    writer.writerow(
        ["run", "iteration", "time", "mean_velocity", "mean_density", "mean_entropy", "mean_internal_energy",
         "mean_pressure", "mean_temperature"])
    for run in runs:
        path = base_path + f"/{run}/formatted_{run}"
        # loop through all files in path, where the iteration is the file name minus the extension
        for file in os.listdir(path):
            f = path + "/" + file
            iteration = file.split(".")[0]
            if int(iteration) >= min_iteration and int(iteration) <= max_iteration:
                time = get_time(f)
                df = pd.read_csv(f, skiprows=2, sep=",")
                # use vx, vy, and vz to calculate the magnitude of the velocity
                df["velocity"] = (df["vx"] ** 2 + df["vy"] ** 2 + df["vz"] ** 2) ** 0.5
                # get the initial conditions for the hydrodynamics
                disk = df[df['end_label'] == "DISK"]
                mean_vel = disk['velocity'].mean()
                mean_density = disk['density'].mean()
                mean_entropy = disk['entropy'].mean()
                mean_internal_energy = disk['internal_energy'].mean()
                mean_pressure = disk['pressure'].mean()
                mean_temperature = disk['temperature'].mean()
                writer.writerow(
                    [run, iteration, time, mean_vel, mean_density, mean_entropy, mean_internal_energy, mean_pressure,
                     mean_temperature])
    outfile.close()


def consolidate_ic():
    # identify the max pressure of each run and write to a file, which you can use to identify the time of impact
    outfile = open(to_dir + "/time_of_impact.csv", 'w')
    writer = csv.writer(outfile, delimiter=",")
    writer.writerow(
        ["run", "iteration", "time", "mean_velocity", "mean_density", "mean_entropy", "mean_internal_energy",
         "mean_pressure", "mean_temperature"])
    for run in runs:
        df = pd.read_csv(to_dir + "/" + run, sep=",")
        max_pressure = df['mean_pressure'].max()
        time_of_impact = df[df['mean_pressure'] == max_pressure]['time'].values[0]
        # get the time of impact row as a list
        row = df[df['time'] == time_of_impact].values[0].tolist()
        # write row to file
        writer.writerow(row)

        # plot the time of impact
        # get the iteration of the time of impact
        iteration = df[df['time'] == time_of_impact]['iteration'].values[0]
        # get the file name of the time of impact
        file = base_path + f"/{run}/formatted_{run}/" + f"{iteration}.csv"
        # read the file
        df = pd.read_csv(file, skiprows=2, sep=",")
        planet, disk, escape = df[df['end_label'] == "PLANET"], df[df['end_label'] == "DISK"], df[
            df['end_label'] == "ESCAPE"]
        # use dark background
        plt.style.use('dark_background')
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(planet['x'], planet['y'], s=1, label="planet")
        ax.scatter(disk['x'], disk['y'], s=1, label="disk")
        ax.scatter(escape['x'], escape['y'], s=1, label="escape")
        # title should include run, time, and iteration
        ax.set_title(f"{run} - time: {time_of_impact} hours, iteration: {iteration}")
        # save the figure as the run name
        plt.savefig(fig_dir + "/" + run + ".png", dpi=200)
    outfile.close()

def plot_all_variables():
    df = pd.read_csv(to_dir + "/disk_history.csv", sep=",").sort_values(by=['iteration'])
    # need a plot with 2 rows and 3 columns
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    # want to plot mean velocity, mean density, mean entropy, mean internal energy, mean pressure, and mean temperature
    # against time
    axs[0, 0].plot(df['time'], df['mean_velocity'])
    axs[0, 0].set_title("mean velocity")
    axs[0, 1].plot(df['time'], df['mean_density'])
    axs[0, 1].set_title("mean density")
    axs[0, 2].plot(df['time'], df['mean_entropy'])
    axs[0, 2].set_title("mean entropy")
    axs[1, 0].plot(df['time'], df['mean_internal_energy'])
    axs[1, 0].set_title("mean internal energy")
    axs[1, 1].plot(df['time'], df['mean_pressure'])
    axs[1, 1].set_title("mean pressure")
    axs[1, 2].plot(df['time'], df['mean_temperature'])
    axs[1, 2].set_title("mean temperature")
    plt.savefig(fig_dir + "/all_variables.png", dpi=200)

get_ic()
plot_all_variables()
consolidate_ic()
