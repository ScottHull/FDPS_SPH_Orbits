import csv
import numpy as np

from new_and_old_eos import vmf_evolution
from src.new_and_old_eos import get_particles
from src.time import get_max_time, seconds_to_hours

min_iteration = 0
max_iteration = 2650
increment = 10
num_processes_new = 200
num_processes_old = 100
new_eos_silicate_phase_curve = "src/phase_data/forstSTS__vapour_curve.txt"
old_eos_silicate_phase_curve = "src/phase_data/duniteN_vapour_curve.txt"
new_eos_formatted_path = "/home/theia/scotthull/1M/formatted_gi_new_eos"
old_eos_formatted_path = "/home/theia/scotthull/1M/formatted_gi_old_eos"
new_eos_unformatted_path = "/home/theia/scotthull/1M/gi_new_eos"
old_eos_unformatted_path = "/home/theia/scotthull/1M/gi_new_eos"

# get the max simulated time
max_time = get_max_time(max_iteration=max_iteration, path=new_eos_formatted_path)

# instantiate profiling classes
vmf = vmf_evolution.VMFtimeseries(
    new_phase_path=new_eos_silicate_phase_curve,
    old_phase_path=old_eos_silicate_phase_curve,
    min_time=0,
    max_time=max_time,
    txt_fname="vmf_b_073.txt",
    output="/home/theia/scotthull/1M/vmf_073",
    animation_name="vmf_timeseries_b_073.mp4"
)

# cycle through outputs
for iteration in np.arange(min_iteration, max_iteration + increment, increment):
    if iteration < 50 and iteration % 10 == 0:
        new_particles, new_time = get_particles(path=new_eos_formatted_path, number_processes=num_processes_new,
                                                time=iteration,
                                                solve=True)
        old_particles, old_time = get_particles(path=old_eos_formatted_path, number_processes=num_processes_old,
                                                time=iteration,
                                                solve=True)

