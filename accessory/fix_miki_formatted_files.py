#!/usr/bin/env python3
import os
import csv
import pandas as pd

def get_time_and_particle_count(f):
    formatted_time = None
    particle_count = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
        particle_count = int(next(reader)[0])
    infile.close()
    return formatted_time, particle_count

to_path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new_high/miki_combined_format"

for file in os.listdir(to_path):
    df = pd.read_csv(to_path + "/" + file, skiprows=2, delimiter=" ")
    for header in df.keys():
        if header not in ["ID", "TYPE"]:
            df[header] = ["{:.8E}".format(i) for i in df[header]]
    fname = str(int(file.replace(".dat", ""))).zfill(4)
    df.to_csv(to_path + "/{}.dat".format(fname), sep=" ", index=False)
    os.remove(to_path + "/" + file)
