#!/usr/bin/env python3
import os
import csv
import numpy as np
import pandas as pd

min_iteration = 0
max_iteration = 3000
number_processes = 200
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
cutoff_densities = [5, 500, 1000, 2000]
runs = ['new', 'old']
angles = ['b073', 'b075']
ftemplate = "results.{}_{}_{}.dat"

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

def get_all_sims(high=False):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for r in runs:
        for a in angles:
            n = "n"
            if r == "old":
                n = "o"
            for cd in cutoff_densities:
                output_name = fformat.format(cd, a, r)
                title_name = tformat.format(cd, a, n)
                titles.append(title_name)
                names.append(output_name)
                if cd == 5 and high and r == "new":
                    output_name = fformat.format(cd, a, r) + "_high"
                    names.append(output_name)
                    title_name = tformat.format(cd, a, n) + "-high"
                    titles.append(title_name)
        return names, titles



sims, titles = get_all_sims(high=False)
data = {}
for sim in sims:
    data.update({sim: []})
    for iteration in np.arange(min_iteration, max_iteration + 1, 1):
        path = base_path + "{}/{}/".format(sim, sim)
        file = path + ftemplate.format(str(iteration).zfill(5),
                                       str(number_processes).zfill(5),
                                       str(0).zfill(5))
        data[sim].append(get_time(file))

df = pd.DataFrame(data, index=list(np.arange(min_iteration, max_iteration + 1, 1)))
df['iteration'] = list(np.arange(min_iteration, max_iteration + 1, 1))
df.to_csv("iterations_vs_times.csv")
