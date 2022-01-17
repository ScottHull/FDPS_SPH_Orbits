#!/usr/bin/env python3
import os
import numpy as np
import csv
import matplotlib.pyplot as plt
from src.report import BuildReports
import multiprocessing as mp

"""
This code produces cumulative outputs, but only "solves" for the disk at certain time intervals.
Useful for getting a coarse understanding of disk mass change without having to solve at each timestep.
"""

plt.style.use("dark_background")

path = "/home/theia/scotthull/1M_high_rho_cutoff/gi_new_eos_b_073_high_rho_cutoff_1M"
to_path = "/home/theia/scotthull/1M_high_rho_cutoff/formatted_gi_new_eos_b_073_high_rho_cutoff_1M"
if "new" in path:
    phase_curve = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/forstSTS__vapour_curve.txt"
else:
    phase_curve = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/duniteN__vapour_curve.txt"
num_processes = 200
assess_at_interval = 100
start_iteration = 0
end_iteration = 3000
increment = 1
end_disk_state_iteration = 3000

os.mkdir(to_path)

reports = BuildReports(
    accessory_path="/home/theia/scotthull/1M_high_rho_cutoff/gi_new_eos_b_073_high_rho_cutoff_1M_disk_reports",
    start_time=start_iteration,
    end_time=end_disk_state_iteration,
    eos_phase_path=phase_curve,
    from_dir=path,
    number_processes=num_processes,
    to_dir=to_path
)

def __loop(time):
    print("At iteration {}".format(time))
    if time % assess_at_interval == 0 and time != 0:
        reports.build_report_at_time(time=time, solve=True)
    else:
        reports.build_report_at_time(time=time, solve=False)

pool = mp.Pool(5)
pool.map(__loop, np.arange(start_iteration, end_iteration + increment, increment))
