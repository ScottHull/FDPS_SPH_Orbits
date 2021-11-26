#!/usr/bin/env python3
import csv
import numpy as np
import matplotlib.pyplot as plt

from new_and_old_eos import vmf_evolution, disk_properties
from src.new_and_old_eos import get_particles
from src.time import get_max_time, seconds_to_hours
from src.read import get_particles_from_formatted

min_iteration = 0
max_iteration = 3000
increment = 20
num_processes_new = 200
num_processes_old = 200
new_eos_silicate_phase_curve = "src/phase_data/forstSTS__vapour_curve.txt"
old_eos_silicate_phase_curve = "src/phase_data/duniteN_vapour_curve.txt"
new_eos_formatted_path = "/home/theia/scotthull/1M/formatted_gi_new_eos_b_073"
old_eos_formatted_path = "/home/theia/scotthull/1M/formatted_gi_old_eos_b_073"
new_eos_unformatted_path = "/home/theia/scotthull/1M/gi_new_eos_b_073"
old_eos_unformatted_path = "/home/theia/scotthull/1M/gi_old_eos_b_073"


# # profile the end-state disk
# parameters = ['entropy', 'pressure', 'density', 'temperature']
# disk_profile = disk_properties.DiskProperties(
#     new_eos_path=new_eos_unformatted_path,
#     old_eos_path=old_eos_unformatted_path,
#     properties=None,
#     formatted=True,
#     new_num_processes=num_processes_new,
#     old_num_processes=num_processes_old
# )
# plt.style.use("dark_background")
# fig, axs = plt.subplots(len(parameters), 2, figsize=(12, 20), sharex="all",
#                                 gridspec_kw={"hspace": 0.10, "wspace": 0.14})
# fig.patch.set_facecolor('xkcd:black')
# pm_end_new_eos, pm_end_old_eos = disk_profile.get_end_state_disk_particles(end_iteration=max_iteration)
# disk_profile.profile_disk_at_time(
#     fig=fig,
#     axs=axs,
#     iteration=max_iteration,
#     new_eos_particles=pm_end_new_eos,
#     old_eos_particles=pm_end_old_eos,
# )



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
    new_time, new_particles = get_particles(
        number_processes=num_processes_new,
        path=new_eos_unformatted_path,
        time=iteration
    )
    old_time, old_particles = get_particles(
        number_processes=num_processes_old,
        path=old_eos_unformatted_path,
        time=iteration
    )


