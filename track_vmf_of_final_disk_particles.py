#!/usr/bin/env python3
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt

from src.vapor import calc_vapor_mass_fraction_with_circularization_from_formatted

runs = ['500_b073_new']

min_iteration = 0
max_iteration = 1800

base_path = "/home/theia/scotthull/Paper1_SPH/gi"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"

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


times, vmfs = [], []
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
            vmf = calc_vapor_mass_fraction_with_circularization_from_formatted(disk, phase_path=new_phase_path,
                                                                               df_label="end_label")
            times.append(time)
            vmfs.append(vmf)

df = pd.DataFrame({"time": times, "vmf": vmfs})
df.to_csv("vmf_of_final_disk_particles.csv", index=False)
