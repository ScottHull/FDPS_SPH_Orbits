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
phase_curve = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/forstSTS__vapour_curve.txt"
num_processes = 200
assess_at_interval = 100
start_iteration = 0
end_iteration = 2900
end_disk_state_iteration = 1000

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
    if time % assess_at_interval:
        reports.build_report_at_time(time=time, solve=True)
    else:
        reports.build_report_at_time(time=time, solve=False)

pool = mp.Pool(5)
for time in np.arange(0, end_iteration + 1, 1):
    pool.map(__loop, [time])
pool.close()
pool.join()


