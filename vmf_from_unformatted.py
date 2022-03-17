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

from src.vapor import calc_vapor_mass_fraction_with_circularization, calc_vapor_mass_fraction
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
    iteration, cd, output_name, path, to_path, to_path2 = arg
    f1 = to_path + "/{}.csv".format(iteration)
    f2 = to_path2 + "/vmf_with_circ_{}_{}_{}.csv".format(angle, runs, iteration)
    if os.path.exists(f2):
        os.remove(f2)
    outfile = open(f2, 'w')
    header = "run,iteration,time,entropy_no_circ_disk,delta_s_circ_disk,total_new_entropy_disk,vmf_no_circ,vmf_circ\n"
    outfile.write(header)
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
    write_report_at_time(particles=particles, fname=f1)

    try:
        mean_s_no_circ = mean([p.entropy for p in particles if p.label == "DISK"])
    except:
        mean_s_no_circ = 0.0
    try:
        mean_delta_s_circ = mean([p.circularization_entropy_delta for p in particles if p.label == "DISK"])
    except:
        mean_delta_s_circ = 0.0
    try:
        mean_total_s = mean([p.entropy + p.circularization_entropy_delta for p in particles if p.label == "DISK"])
    except:
        mean_total_s = 0.0
    vmf_no_circ = calc_vapor_mass_fraction(particles=particles, phase_path=phase_path)
    vmf_circ = calc_vapor_mass_fraction_with_circularization(particles=particles, phase_path=phase_path)
    line = "{},{},{},{},{},{},{},{}\n".format(
        output_name, iteration, formatted_time, mean_s_no_circ, mean_delta_s_circ, mean_total_s, vmf_no_circ, vmf_circ
    )
    outfile.write(line)
    outfile.close()
    return 0

def build_reports():
    for cd in cutoff_densities:
        output_name = "{}_{}_{}".format(cd, angle, runs)
        path = base_path + "{}_{}_{}/{}_{}_{}".format(cd, angle, runs, cd, angle, runs)
        to_path = base_path + output_name + "/circularized_{}".format(output_name)
        to_path2 = base_path + output_name + "/circularized_{}_disk_descriptions".format(output_name)
        if not os.path.exists(to_path):
            os.mkdir(to_path)
        if not os.path.exists(to_path2):
            os.mkdir(to_path2)
        pool = mp.Pool(5)
        pool.map(mp_task, [[iteration, cd, output_name, path, to_path, to_path2] for iteration in
                           list(np.arange(min_iteration, max_iteration + increment, increment))])
        pool.close()
        pool.join()

def add_times():
    for runs in ["new", "old"]:
        for cd in cutoff_densities:
            for iteration in np.arange(min_iteration, max_iteration + increment, increment):
                print("at {} // {} // {}".format(runs, cd, iteration))
                path = base_path + "{}_{}_{}/{}_{}_{}".format(cd, angle, runs, cd, angle, runs)
                to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
                cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
                combined_file = cf.combine()
                formatted_time = round(cf.sim_time * 0.000277778, 2)
                os.remove(to_fname)
                output_name = "{}_{}_{}".format(cd, angle, runs)
                to_path2 = base_path + output_name + "/circularized_{}_disk_descriptions".format(output_name)
                f2 = to_path2 + "/vmf_with_circ_{}_{}_{}.csv".format(angle, runs, iteration)
                df = pd.read_csv(f2)
                df['time'] = [formatted_time]
                df.to_csv(f2)


def plot_vmfs():
    fig, axs = plt.subplots(2, 2, figsize=(16, 9), sharex='all', sharey='all')
    axs = axs.flatten()
    for runs in ["new", "old"]:
        r = "n"
        vmf_total_index = 0
        vmf_no_circ_index = 1
        if runs == "old":
            r = "o"
            vmf_total_index += 2
            vmf_no_circ_index += 2
        for cd in cutoff_densities:
            times = []
            vmfs_total = []
            vmfs_without_circ = []
            delta_s_circ = []
            total_ss = []
            for iteration in np.arange(min_iteration, max_iteration + increment, increment):
                print("at {} // {} // {}".format(runs, cd, iteration))
                output_name = "{}_{}_{}".format(cd, angle, runs)
                to_path2 = base_path + output_name + "/circularized_{}_disk_descriptions".format(output_name)
                f2 = to_path2 + "/vmf_with_circ_{}_{}_{}.csv".format(angle, runs, iteration)
                df = pd.read_csv(f2)
                formatted_time = df['time']
                vmf_no_circ = df['vmf_no_circ']
                vmf_total = df['vmf_circ']
                total_s = df['total_new_entropy_disk']
                delta_s_due_to_circ = df['delta_s_circ_disk']
                s_no_circ = df['entropy_no_circ_disk']
                total_ss.append(total_s)
                delta_s_circ.append(delta_s_due_to_circ)

                times.append(formatted_time)
                vmfs_total.append(vmf_total * 100)
                vmfs_without_circ.append(vmf_no_circ * 100)
            axs[vmf_total_index].plot(
                times, vmfs_total, linewidth=2.0, label="{}{}{}".format(cd, angle, runs)
            )
            axs[vmf_no_circ_index].plot(
                times, vmfs_without_circ, linewidth=2.0, label="{}{}{}".format(cd, angle, runs)
            )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.legend(loc='upper right')
    axs[0].set_title("VMF (with orbital circularization)")
    axs[1].set_title("VMF (without orbital circularization)")
    axs[0].set_ylabel("VMF (%)")
    axs[2].set_ylabel("VMF (%)")
    axs[2].set_xlabel("Time (hrs)")
    axs[3].set_xlabel("Time (hrs)")
    plt.savefig("vmf_w_wo_circ.png", format='png', dpi=200)

    fig, axs = plt.subplots(2, 2, figsize=(16, 9), sharex='all')
    axs = axs.flatten()
    for runs in ["new", "old"]:
        r = "n"
        vmf_total_index = 0
        vmf_no_circ_index = 1
        if runs == "old":
            r = "o"
            vmf_total_index += 2
            vmf_no_circ_index += 2
        for cd in cutoff_densities:
            times = []
            vmfs_total = []
            vmfs_without_circ = []
            delta_s_circ = []
            total_ss = []
            for iteration in np.arange(min_iteration, max_iteration + increment, increment):
                print("at {} // {} // {}".format(runs, cd, iteration))
                output_name = "{}_{}_{}".format(cd, angle, runs)
                to_path2 = base_path + output_name + "/circularized_{}_disk_descriptions".format(output_name)
                f2 = to_path2 + "/vmf_with_circ_{}_{}_{}.csv".format(angle, runs, iteration)
                df = pd.read_csv(f2)
                formatted_time = df['time']
                vmf_no_circ = df['vmf_no_circ']
                vmf_total = df['vmf_circ']
                total_s = df['total_new_entropy_disk']
                delta_s_due_to_circ = df['delta_s_circ_disk']
                s_no_circ = df['entropy_no_circ_disk']
                total_ss.append(total_s)
                delta_s_circ.append(delta_s_due_to_circ)

                times.append(formatted_time)
                vmfs_total.append(vmf_total * 100)
                vmfs_without_circ.append(vmf_no_circ * 100)
            axs[vmf_total_index].plot(
                times, total_ss, linewidth=2.0, label="{}{}{}".format(cd, angle, runs)
            )
            axs[vmf_no_circ_index].plot(
                times, delta_s_circ, linewidth=2.0, label="{}{}{}".format(cd, angle, runs)
            )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.legend(loc='upper right')
    axs[0].set_title("Avg. Entropy (with orbital circularization)")
    axs[1].set_title(r"Avg. $\Delta$S Due to Orbital Circularization")
    axs[0].set_ylabel("S")
    axs[2].set_ylabel("S")
    axs[2].set_xlabel("Time (hrs)")
    axs[3].set_xlabel("Time (hrs)")
    plt.savefig("s_w_wo_circ.png", format='png', dpi=200)


def fix_entropies():
    for runs in ["new", "old"]:
        for cd in cutoff_densities:
            for iteration in np.arange(min_iteration, max_iteration + increment, increment):
                print("at {} // {} // {}".format(runs, cd, iteration))
                output_name = "{}_{}_{}".format(cd, angle, runs)
                to_path = base_path + output_name + "/circularized_{}".format(output_name)
                to_path2 = base_path + output_name + "/circularized_{}_disk_descriptions".format(output_name)
                f1 = to_path + "/{}.csv".format(iteration)
                f2 = to_path2 + "/vmf_with_circ_{}_{}_{}.csv".format(angle, runs, iteration)
                df = pd.read_csv(f1)
                disk = df[df['label'] == "DISK"]
                df2 = pd.read_csv(f2)
                mean_delta_s = mean(disk['circ_entropy_delta'])
                mean_total_s_no_circ = mean(disk['entropy'])
                mean_total_s_w_circ = mean(disk['circ_entropy_delta'] + disk['entropy'])

                df2['total_new_entropy_disk'] = mean_total_s_w_circ
                df2['delta_s_circ_disk'] = mean_delta_s
                df2['entropy_no_circ_disk'] = mean_total_s_no_circ
                df2.to_csv(f2)



def fix_delimiting():
    def replace_str_index(text, index=0, replacement=''):
        return '%s%s%s' % (text[:index], replacement, text[index + 1:])

    find_all = lambda c, s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    for runs in ["new", "old"]:
        for cd in cutoff_densities:
            for iteration in np.arange(min_iteration, max_iteration + increment, increment):
                print("at {} // {} // {}".format(runs, cd, iteration))
                path = base_path + "{}_{}_{}/{}_{}_{}".format(cd, angle, runs, cd, angle, runs)
                output_name = "{}_{}_{}".format(cd, angle, runs)
                to_path2 = base_path + output_name + "/circularized_{}_disk_descriptions".format(output_name)
                f2 = to_path2 + "/vmf_with_circ_{}_{}_{}.csv".format(angle, runs, iteration)
                tmp_p = to_path2 + "/tmp.csv"
                if os.path.exists(tmp_p):
                    os.remove(tmp_p)
                with open(tmp_p, 'w') as outfile:
                    with open(f2, 'r') as infile:
                        header = next(infile)
                        outfile.write(header)
                        for row in infile:
                            periods = find_all(str(row), ".")
                            row = replace_str_index(str(row), periods[-3], ',')
                            outfile.write(row)
                    infile.close()
                outfile.close()
                os.remove(f2)
                os.rename(tmp_p, f2)


# build_reports()
plot_vmfs()


