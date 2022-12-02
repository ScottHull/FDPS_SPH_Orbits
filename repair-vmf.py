#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd

from src.interpolation import NearestNeighbor1D

runs = [
    '5_b073_new', '5_b073_old', '5_b073_new_high', '500_b073_new', '500_b073_old', '1000_b073_new', '1000_b073_old',
    '2000_b073_new', '2000_b073_old',
    '5_b075_new', '5_b075_old', '500_b075_new', '500_b075_old',
    '1000_b075_new', '1000_b075_old', '2000_b075_new', '2000_b075_old', '2000_b075_old_low'
]
min_iteration = 0
max_iteration = 1800
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                       names=["temperature", "density_sol_liq", "density_vap", "pressure",
                              "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                       names=["temperature", "density_sol_liq", "density_vap", "pressure",
                              "entropy_sol_liq", "entropy_vap"])

def calc_vmf(data_df, phase_df, entropy_header):
    disk_df = data_df[data_df['label'] == "DISK"]
    disk_df = disk_df[disk_df['tag'] % 2 == 0]
    disk_df['entropy_w_circ'] = disk_df['entropy'] + disk_df['circ_entropy_delta']
    disk_df = disk_df[disk_df['circ_entropy_delta'] < 5000]

    nearest_neighbor = NearestNeighbor1D()
    supercritical = max(phase_df['temperature'])
    # get the location of the phase curve that corresponds to the supercritical point
    supercritical_idx = phase_df[phase_df['temperature'] == supercritical].index[0]
    # get the entropy values for the supercritical point
    supercritical_entropy = phase_df.iloc[supercritical_idx]['entropy_sol_liq']
    vmf = 0.0
    num_particles = 0
    # particles = {'liq': [], 'vap': [], 'mix': [], 'supercritical': []}
    for s, t in zip(disk_df[entropy_header], data_df['temperature']):
        nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                    array=list(phase_df['temperature']))
        entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
        entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
        if t >= supercritical:
            vmf += 1.0
            # particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
        elif s < entropy_liq:
            vmf += 0.0
            # particles['liq'].append((s, t, entropy_liq, entropy_vap))
        elif entropy_liq <= s <= entropy_vap:
            vmf += (s - entropy_liq) / (entropy_vap - entropy_liq)
            # particles['mix'].append((s, t, entropy_liq, entropy_vap))
        elif s > entropy_vap:
            vmf += 1.0
            # particles['vap'].append((s, t, entropy_liq, entropy_vap))
        num_particles += 1
    try:
        return vmf / num_particles
    except ZeroDivisionError:
        return 0.0

for run in runs:
    increment = 50
    if "high" in run:
        increment = 100
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        try:
            phase_df = new_phase_df if 'new' in run else old_phase_df
            circ_path = base_path + f"{run}/circularized_{run}/{iteration}.csv"
            report_path = base_path + f"{run}/{run}_reports/{iteration}.csv"
            report_path_2 = base_path + f"{run}/{run}_reports2"
            if not os.path.exists(report_path_2):
                os.mkdir(report_path_2)
            circ_df = pd.read_csv(circ_path)
            circ_df['entropy_w_circ'] = circ_df['entropy'] + circ_df['circ_entropy_delta']
            vmf_w_circ = calc_vmf(circ_df, phase_df, 'entropy_w_circ')
            vmf_wo_circ = calc_vmf(circ_df, phase_df, 'entropy')
            report_df = pd.read_csv(report_path)
            report_df['DISK_VMF_W_CIRC'] = vmf_w_circ
            report_df['DISK_VMF_WITHOUT_CIRC'] = vmf_wo_circ
            report_df.to_csv(report_path_2 + f"/{iteration}.csv", index=False)
        except:
            pass


