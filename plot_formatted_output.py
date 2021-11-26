#!/usr/bin/env python3
import pandas
import matplotlib.pyplot as plt

from src.new_and_old_eos import get_particles, get_parameter_from_particles

iteration = 3000
num_processes_new = 200
num_processes_old = 200
new_eos_formatted_path = "/home/theia/scotthull/1M/formatted_gi_new_eos_b_073"
old_eos_formatted_path = "/home/theia/scotthull/1M/formatted_gi_old_eos_b_073"
r_earth = 6371 * 1000

labels = {
    1: "Target Silicate",
    2: "Target Iron",
    3: "Impactor Silicate",
    4: "Impactor Iron"
}

plt.style.use("dark_background")

vars = ["entropy", "temperature", "pressure"]

fig, axs = plt.subplots(len(vars), 2, figsize=(8, 16), sharex='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.05})
fig.patch.set_facecolor('xkcd:black')

new_particles, new_time = get_particles(
    number_processes=num_processes_new,
    path=new_eos_formatted_path,
    time=iteration,
    solve=False,  # if using formatted, already solved!
    formatted=True
)
old_particles, old_time = get_particles(
    number_processes=num_processes_old,
    path=old_eos_formatted_path,
    time=iteration,
    solve=False,  # if using formatted, already solved!
    formatted=True
)

plotting_index = 0
for v in vars:
    new_index, old_index = plotting_index, plotting_index + 1
    new_ax, old_ax = axs.flatten()[new_index], axs.flatten()[old_index]
    both = [new_ax, old_ax]
    all_y_values = []
    for tag in labels.keys():
        new_y = [get_parameter_from_particles(particle=p, parameter=v) for p in new_particles if p.tag == tag]
        old_y = [get_parameter_from_particles(particle=p, parameter=v) for p in old_particles if p.tag == tag]
        all_y_values += new_y + old_y
        new_ax.scatter(
            [p.distance / r_earth for p in new_particles if p.tag == tag],
            new_y,
            s=2,
            label=labels[tag],
        )
        old_ax.scatter(
            [p.distance / r_earth for p in old_particles if p.tag == tag],
            old_y,
            s=2,
            label=labels[tag]
        )
    old_ax.set_yticks([])
    old_ax.set_yticks([], minor=True)
    new_ax.set_ylabel(v)
    for ax in both:
        ax.set_ylim(min(all_y_values), max(all_y_values))
        ax.grid(alpha=0.4)
    if v == vars[-1]:
        for ax in both:
            ax.set_xlabel(r'Radius $R_{\bigoplus}$')
    plotting_index += 2

axs.flatten()[0].set_title("New EoS")
axs.flatten()[1].set_title("Old EoS")
legend = axs.flatten()[0].legend(loc='upper right')
for handle in legend.legendHandles:
    handle.set_sizes([4.0])
plt.savefig("end_state_formatted.png", format='png')
