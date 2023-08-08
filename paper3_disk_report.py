#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18,})
plt.style.use("dark_background")

from src.combine import CombineFile
from src.identify import ParticleMap
from src.report import get_sim_report, write_report_at_time, build_latex_table_from_disk_report, rows_map

base_path = "/home/theia/scotthull/Paper2_SPH/gi/"
sim = ["500_mars"]
min_iteration = 0
max_iteration = 1800
increment = 50
number_processes = 200

# MASS_EARTH = 5.972e24  # kg
# RADIUS_EARTH = 6.371e6  # m
MASS_MARS = 6.39e23  # kg
RADIUS_MARS = 3.389e6  # m

# CHANGE THESE
MASS_PROTOPLANET = MASS_MARS
RADIUS_PROTOPLANET = RADIUS_MARS

phase_path = "src/phase_data/forstSTS__vapour_curve.txt"

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

def get_all_sims(**kwargs):
    return [sim[0]], [sim[1]]


def __build_report(args):
    iteration, output_name, title = args
    to_path = base_path + output_name + "/circularized_{}".format(output_name)
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    f1 = to_path + "/{}.csv".format(iteration)
    path = base_path + "{}/{}".format(output_name, output_name)
    to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
    combined_file = cf.combine()
    formatted_time = round(cf.sim_time * 0.000277778, 2)
    f = os.getcwd() + "/{}".format(to_fname)
    pm = ParticleMap(path=f, center=True, relative_velocity=False, mass_protoplanet=MASS_PROTOPLANET,
                     initial_equatorial_radius=RADIUS_PROTOPLANET + 1000)
    particles = pm.collect_particles(find_orbital_elements=True)
    os.remove(to_fname)
    pm.solve(particles=particles, phase_path=phase_path, report_name="{}-report.txt".format(output_name),
             iteration=iteration, simulation_time=formatted_time)
    to_report_path = base_path + "{}/{}_reports".format(output_name, output_name)
    write_report_at_time(particles=particles, fname=f1)
    get_sim_report(particle_df=pd.read_csv(f1), phase_path=phase_path, formatted_time=formatted_time,
                   iteration=iteration, sim_name=title, to_path=to_report_path)

def __base_build_report(args):
    to_path, iteration, output_name, titles, index = args
    f1 = to_path + "/{}.csv".format(iteration)
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


def __build_report_iteration_mp(args):
    index, output_name, titles = args
    to_path = base_path + output_name + "/circularized_{}".format(output_name)
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    pool = mp.Pool(5)
    pool.map(__base_build_report, [[to_path, iteration, output_name, titles, index] for iteration in np.arange(min_iteration,
                                                                          max_iteration + increment, increment)])
    pool.close()
    pool.join()


def __plot_disk_report(run_names: list, run_titles: list, to_base_path: str, imp_angle: str, iteration: int):
    def __map_run_name(r):
        if "5_" in r.replace("b075", ""):
            return 0
        if "500_" in r.replace("b075", ""):
            return 1
        if "1000_" in r.replace("b075", ""):
            return 2
        if "2000_" in r.replace("b075", ""):
            return 3

    headers = ["MEAN_DISK_ENTROPY_W_CIRC", "DISK_VMF_W_CIRC", "DISK_MASS", "DISK_ANGULAR_MOMENTUM", "DISK_THEIA_MASS_FRACTION",
               "DISK_IRON_MASS_FRACTION"]
    cutoff_densities = [5, 500, 1000, 2000]
    fig, axs = plt.subplots(int(len(headers) / 2), 2, figsize=(16, 24), sharex="all",
                            gridspec_kw={"hspace": 0.14, "wspace": 0.14})
    axs = axs.flatten()

    for h_index, h in enumerate(headers):
        y_new, y_old = [None, None, None, None], [None, None, None, None]
        for index, run in enumerate(run_names):
            path = to_base_path + "{}/{}_reports/".format(run, run)
            df = pd.read_csv(path + "{}.csv".format(iteration))
            val = float(str(df[h][0]).split(" ")[0])
            if "new" in run:
                y_new[__map_run_name(run)] = val
            else:
                y_old[__map_run_name(run)] = val
        axs[h_index].plot(
            cutoff_densities, y_new, linewidth=2.0, label=r"Stewart M-ANEOS ($N = {10^6}$)"
        )
        axs[h_index].plot(
            cutoff_densities, y_old, linewidth=2.0, label=r"N-SPH M-ANEOS ($N = {10^6}$)"
        )
        axs[-2].set_xlabel("Cutoff Density")
        axs[-1].set_xlabel("Cutoff Density")
        label = rows_map[h][1:-1]
        axs[h_index].set_ylabel(label)
        axs[h_index].grid(alpha=0.4)
        axs[h_index].legend(loc='upper left')
    plt.savefig("{}_disk_report.png".format(imp_angle), format='png', dpi=200)



def build_report():
    sims, titles = get_all_sims(high=False)
    pool = mp.Pool(5)
    pool.map(__build_report, [[i, sims[0], titles[0]] for i in np.arange(min_iteration, max_iteration + increment, increment)])
    pool.close()
    pool.join()

def build_report_high_res():
    sims, titles = ['5_b073_new_high'], ['5b073S-high']
    __build_report_iteration_mp(args=[0, sims[0], titles])

def make_report_latex_table():
    sims, titles = get_all_sims(high=False)
    build_latex_table_from_disk_report(run_names=sims, run_titles=titles, to_base_path=base_path,
                            filename="{}_disk_latex_table.txt".format(angle), iteration=max_iteration)

def plot_disk_report():
    sims, titles = get_all_sims(high=False)
    __plot_disk_report(run_names=sims, run_titles=titles, to_base_path=base_path, imp_angle=angle,
                       iteration=max_iteration)
