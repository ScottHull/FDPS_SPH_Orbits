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


times, vmfs = [], []
end_time_df = pd.read_csv(
    os.path.join(base_path, f"{max_iteration}.csv"),
)
end_time_disk = end_time_df[end_time_df["label"] == "DISK"]
end_time_particle_ids = end_time_disk["id"].values
for run in runs:
    path = base_path + f"/{run}/circularized_{run}"
    path2 = base_path + f"/{run}/{run}_reports"
    # loop through all files in path, where the iteration is the file name minus the extension
    for file in os.listdir(path):
        f = path + "/" + file
        iteration = file.split(".")[0]
        if int(iteration) >= min_iteration and int(iteration) <= max_iteration:
            # read the report file as a pandas df and get the time
            f2 = path2 + "/" + iteration + ".csv"
            df2 = pd.read_csv(f2, sep=",")
            time = df2["TIME_HRS"][0]
            df = pd.read_csv(f, sep=",")
            # use vx, vy, and vz to calculate the magnitude of the velocity
            df["velocity"] = (df["vx"] ** 2 + df["vy"] ** 2 + df["vz"] ** 2) ** 0.5
            # get the initial conditions for the hydrodynamics
            # limit df to only the particles that were in the disk at the end of the simulation
            disk = df[df["id"].isin(end_time_particle_ids)]
            vmf = calc_vapor_mass_fraction_with_circularization_from_formatted(disk, phase_path=new_phase_path,
                                                                               df_label="end_label")
            times.append(time)
            vmfs.append(vmf)

df = pd.DataFrame({"time": times, "vmf": vmfs})
df.to_csv("vmf_of_final_disk_particles.csv", index=False)
