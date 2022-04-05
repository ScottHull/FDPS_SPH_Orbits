#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import multiprocessing as mp

from src.combine import CombineFile
from src.identify import ParticleMap
from src.report import get_sim_report, write_report_at_time, build_latex_table_from_disk_report

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5, 500, 1000, 2000]
min_iteration = 0
max_iteration = 1800
increment = 50
number_processes = 200

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

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

def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "n"
        if runs == "old":
            n = "o"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
            if cd == 5 and high and runs == "new":
                output_name = fformat.format(cd, angle, runs) + "_high"
                names.append(output_name)
                title_name = tformat.format(cd, angle, n) + "-high"
                titles.append(title_name)
    return names, titles


def __build_report(args):
    index, output_name, titles = args
    to_path = base_path + output_name + "/circularized_{}".format(output_name)
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        f1 = to_path + "/{}.csv".format(iteration)
        phase_path = new_phase_path
        if "old" in output_name:
            phase_path = old_phase_path
        title = titles[index]
        path = base_path + "{}/{}".format(output_name, output_name)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        f = os.getcwd() + "/{}".format(to_fname)
        pm = ParticleMap(path=f, center=True, relative_velocity=False)
        particles = pm.collect_particles(find_orbital_elements=True)
        os.remove(to_fname)
        pm.solve(particles=particles, phase_path=phase_path, report_name="{}-report.txt".format(output_name),
                 iteration=iteration, simulation_time=formatted_time)
        to_report_path = base_path + "{}/{}_reports".format(output_name, output_name)
        write_report_at_time(particles=particles, fname=f1)
        get_sim_report(particle_df=pd.read_csv(f1), phase_path=phase_path, formatted_time=formatted_time,
                       iteration=iteration, sim_name=title, to_path=to_report_path)

def build_report():
    sims, titles = get_all_sims(high=False)
    pool = mp.Pool(5)
    pool.map(__build_report, [[index, output_name, titles] for index, output_name in enumerate(sims)])
    pool.close()
    pool.join()

def make_report_latex_table():
    sims, titles = get_all_sims(high=False)
    build_latex_table_from_disk_report(run_names=sims, run_titles=titles, to_base_path=base_path,
                            filename="{}_disk_latex_table.txt".format(angle), iteration=max_iteration)