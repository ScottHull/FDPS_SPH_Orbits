#!/usr/bin/env python3
import os
import sys
import pandas as pd
from statistics import mean

from src.vapor import calc_vapor_mass_fraction_from_formatted, \
    calc_vapor_mass_fraction_with_circularization_from_formatted

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angles = ["b073", "b075"]
cutoff_densities = [5, 500, 1000, 2000]

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"


def __mean(l):
    try:
        return mean(l)
    except:
        return 0.0


def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    high_res_name = None
    high_res_title = None
    for angle in angles:
        for runs in ["new", "old"]:
            n = "S"
            if runs == "old":
                n = "N"
            for cd in cutoff_densities:
                output_name = fformat.format(cd, angle, runs)
                title_name = tformat.format(cd, angle, n)
                titles.append(title_name)
                names.append(output_name)
                if cd == 5 and high and runs == "new":
                    high_res_name = fformat.format(cd, angle, runs) + "_high"
                    high_res_title = tformat.format(cd, angle, n) + "-high"
        if high_res_name is not None and high_res_title is not None:
            names.append(high_res_name)
            titles.append(high_res_title)
    return names, titles


def reformat():
    sims, titles = get_all_sims(high=False)
    for sim, title in zip(sims, titles):
        phase_path = new_phase_path
        if "old" in sim:
            phase_path = old_phase_path
        formatted_path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        for i in os.listdir(report_path):
            iteration = i.replace(".csv", "")
            df_formatted = pd.read_csv(formatted_path + i)
            df_report = pd.read_csv(report_path + i)
            try:
                del df_report['Unnamed: 0']
                vmf_uncirc = calc_vapor_mass_fraction_from_formatted(df=df_formatted, phase_path=phase_path)
                vmf_circ = calc_vapor_mass_fraction_with_circularization_from_formatted(df=df_formatted,
                                                                                        phase_path=phase_path)
                del df_report['DISK VMF']
                df_report['DISK_VMF_W_CIRC'] = [str(vmf_circ * 100.0) + " %"]
                df_report['DISK_VMF_WITHOUT_CIRC'] = [str(vmf_uncirc * 100.0) + " %"]
                to_f = report_path + "{}.csv".format(iteration)
                df_report.to_csv(to_f, index=False)
                print("Rewrote report at {}".format(to_f))
            except Exception as e:
                print(e)
                pass

def fix_entropies():
    sims, titles = get_all_sims(high=False)
    for sim, title in zip(sims, titles):
        phase_path = new_phase_path
        if "old" in sim:
            phase_path = old_phase_path
        formatted_path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        for i in os.listdir(report_path):
            iteration = i.replace(".csv", "")
            df_formatted = pd.read_csv(formatted_path + i)
            df_report = pd.read_csv(report_path + i)
            try:
                # del df_report['Unnamed: 0']
                del df_report['MEAN_DISK_ENTROPY']
                disk = df_formatted[df_formatted['label'] == "DISK"]
                disk_filtered = disk[disk['circ_entropy_delta'] < 5000]
                df_report['MEAN_DISK_ENTROPY_W_CIRC'] = [str(__mean(disk_filtered['entropy'] + disk_filtered['circ_entropy_delta']))]
                df_report['MEAN_DISK_ENTROPY_WITHOUT_CIRC'] = [str(__mean(disk['entropy']))]
                to_f = report_path + "{}.csv".format(iteration)
                df_report.to_csv(to_f, index=False)
                print("Rewrote report at {}".format(to_f))
            except Exception as e:
                print(e)
                pass

def __add_mean_disk_temp():
    sims, titles = get_all_sims(high=False)
    for sim, title in zip(sims, titles):
        phase_path = new_phase_path
        if "old" in sim:
            phase_path = old_phase_path
        formatted_path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        for i in os.listdir(report_path):
            iteration = i.replace(".csv", "")
            df_formatted = pd.read_csv(formatted_path + i)
            df_report = pd.read_csv(report_path + i)
            try:
                # del df_report['Unnamed: 0']
                disk = df_formatted[df_formatted['label'] == "DISK"]
                df_report['MEAN_DISK_TEMPERATURE'] = [__mean(disk['temperature'])]
                to_f = report_path + "{}.csv".format(iteration)
                df_report.to_csv(to_f, index=False)
                print("Rewrote report at {}".format(to_f))
            except Exception as e:
                print(e)
                pass

# reformat()
# fix_entropies()
__add_mean_disk_temp()
