#!/usr/bin/env python3
import os
import csv
import shutil
from math import pi, asin, isnan
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp

from src.vapor import calc_vapor_mass_fraction_with_circularization
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_max_time, seconds_to_hours
from src.report import write_report_at_time

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
runs = "new"
angle = "b073"
cutoff_densities = [5, 500, 1000, 2000]
min_iteration = 0
max_iteration = 1500
increment = 100
number_processes = 200

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

phase_path = new_phase_path
if runs == "old":
    phase_path = old_phase_path


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


def mp_task(arg):
    iteration, cd, output_name, path, to_path, outfile = arg
    to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    f = os.getcwd() + "/{}".format(to_fname)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=True)
    os.remove(to_fname)
    pm.solve(particles=particles, phase_path=phase_path, report_name="{}-report.txt".format(output_name),
             iteration=iteration, simulation_time=formatted_time)
    write_report_at_time(particles=particles, fname=to_path + "/{}.csv".format(output_name))
    mean_s_no_circ = mean([p.entropy for p in particles if p.label == "DISK"])
    mean_delta_s_circ = mean([p.circularization_entropy_delta for p in particles if p.label == "DISK"])
    mean_total_s = mean([p.entropy + p.circularization_entropy_delta for p in particles if p.label == "DISK"])
    vmf_no_circ = calc_vapor_mass_fraction_with_circularization(particles=particles, phase_path=phase_path)
    vmf_circ = calc_vapor_mass_fraction_with_circularization(particles=particles, phase_path=phase_path)
    line = "{},{},{},{},{}.{},{}\n".format(
        output_name, iteration, mean_s_no_circ, mean_delta_s_circ, mean_total_s, vmf_no_circ, vmf_circ
    )
    outfile.write(line)
    return 0


outfile = open("vmf_with_circ_{}_{}.csv".format(angle, runs), 'w')
header = "run,iteration,entropy_no_circ_disk,delta_s_circ_disk,total_new_entropy_disk,vmf_no_circ,vmf_circ\n"
outfile.write(header)
for cd in cutoff_densities:
    output_name = "{}_{}_{}".format(cd, angle, runs)
    path = base_path + "{}_{}_{}/{}_{}_{}".format(cd, angle, runs, cd, angle, runs)
    to_path = base_path + output_name + "/circularized_{}".format(output_name)
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    pool = mp.Pool(5)
    pool.map(mp_task, [[iteration, cd, output_name, path, to_path, outfile] for iteration in
                       list(np.arange(min_iteration, max_iteration + increment, increment))])
    pool.close()
    pool.join()

outfile.close()
