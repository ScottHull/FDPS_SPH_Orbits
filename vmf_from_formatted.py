import numpy as np
import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt

from src.time import get_max_time, seconds_to_hours
from src.vapor import calc_vapor_mass_fraction_from_formatted


path = "/home/theia/scotthull/1M_500kgm3_cutoff/formatted_gi_new_eos_b_073_500kgm3_end_2900"
phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
start_time = 0
end_time = 2900
increment = 100

plt.style.use("dark_background")
max_time = seconds_to_hours(get_max_time(max_iteration=end_time, path=path))
times = []
vmfs = []
disk_particle_count = []
avg_disk_entropy = []
for time in np.arange(start_time, end_time + increment, increment):
    f = path + "/{}.csv".format(time)
    times.append(get_max_time(max_iteration=time, path=path))
    df = pd.read_csv(f, skiprows=2)
    vmf = calc_vapor_mass_fraction_from_formatted(df=df, phase_path=phase_path) * 100.0
    vmfs.append(vmf)

    disk_particles = df[df['label'] == "DISK"]
    try:
        avg_disk_entropy_at_time = mean(disk_particles['entropy'])
    except:
        avg_disk_entropy_at_time = 0
    avg_disk_entropy.append(avg_disk_entropy_at_time)
    num_disk_particles = len(disk_particles['entropy'])
    disk_particle_count.append(num_disk_particles)
    print(
        "TIME: {}\nVMF: {}\nAVG DISK ENTROPY: {}\nDISK PARTICLE COUNT: {}".format(time, vmf, avg_disk_entropy_at_time,
                                                                                  num_disk_particles)
    )



fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharex="all", gridspec_kw={"wspace": 0.14})
axs = axs.flatten()
ax1, ax2, ax3 = axs
for ax in axs:
    ax.set_xlim(0, max_time)
    ax.grid(alpha=0.4)
    ax.set_xlabel("Time (hrs)")
ax1.plot(
    times, vmfs, linewidth=2.0
)
ax2.plot(
    times, avg_disk_entropy, linewidth=2.0
)
ax3.plot(
    times, disk_particle_count, linewidth=2.0
)
ax1.set_ylabel("VMF (%)"), ax2.set_ylabel("Avg. Disk Entropy"), ax3.set_ylabel("# Disk Particles")
plt.savefig("vmf_from_formatted.png", format='png')
