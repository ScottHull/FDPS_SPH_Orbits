#!/usr/bin/env python3
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt

from src.new_and_old_eos import get_particles

"""
We want to model the iron layer after the GI as a spherical shell.
Volume of a sphere: 4/3 pi r^3
Volume of a spherical shell: 4/3 pi (r_2 - r_1)^3
We need to calculate r_2 - r_1 (can't get it directly from simulation)
Know that density is rho = M / V, so the density of the shell is rho = M / [4/3 pi (r_2 - r_1)^3]
We can get M from the simulation, but what about density?

THERMAL EXPANSIVITY (alpha) : [1 / K]
TEMPERATURE (T) : [K]
DENSITY (rho) : [kg/m3]

DELTA RHO = RHO_0 * ALPHA * DELTA T
RHO_0_iron = 7500

"""

time = 3000
start_time = 0
new_path = "/home/theia/scotthull/1M/formatted_gi_new_eos_b_073"
old_path = "/home/theia/scotthull/1M/formatted_gi_old_eos_b_073"
new_path_unformatted = "/home/theia/scotthull/1M/gi_new_eos_b_073"
old_path_unformatted = "/home/theia/scotthull/1M/gi_old_eos_b_073"
new_processes = 200
old_processes = 200

plt.style.use("dark_background")

labels = {
    0: "Target Silicate",
    1: "Target Iron",
    2: "Impactor Silicate",
    3: "Impactor Iron"
}


def map_t0_particles(t_0_particles):
    mapped = {}
    for p in t_0_particles:
        mapped.update({p.particle_id: p})
    return mapped

def iron_layer_delta_temperature(t_particles, t_0_particles):
    t_0_particles = map_t0_particles(t_0_particles=t_0_particles)
    delta_t = [p.temperature - t_0_particles[p.particle_id].temperature for p in t_particles if p.distance / (6371 * 1000) < 1.5 and p.tag == 3]
    return delta_t

def iron_layer_density_delta(t_particles, t_0_particles):
    t_0_particles = map_t0_particles(t_0_particles=t_0_particles)
    delta_rho = [p.density - t_0_particles[p.particle_id].density for p in t_particles if p.distance / (6371 * 1000) < 1.5 and p.tag == 3]
    return delta_rho


def identify_iron_layer_particles(new_particles, old_particles, earth_radius=6371 * 1000):
    new_impactor_iron_particles = [p for p in new_particles if p.label == "PLANET" and p.tag == 3 and p.distance / (6371 * 1000) < 1.5]
    old_impactor_iron_particles = [p for p in old_particles if p.label == "PLANET" and p.tag == 3 and p.distance / (6371 * 1000) < 1.5]

    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.0, "wspace": 0.10})
    fig.patch.set_facecolor('xkcd:black')
    new_ax, old_ax = axs.flatten()[0], axs.flatten()[1]
    new_ax.scatter(
        [p.distance / earth_radius for p in new_particles if p.distance / (6371 * 1000) < 1.5],
        [p.temperature for p in new_particles if p.distance / (6371 * 1000) < 1.5],
        s=2,
        color='magenta'
    )
    old_ax.scatter(
        [p.distance / earth_radius for p in new_impactor_iron_particles if p.distance / (6371 * 1000) < 1.5],
        [p.temperature for p in old_impactor_iron_particles if p.distance / (6371 * 1000) < 1.5],
        s=2,
        color='magenta'
    )
    new_ax.set_ylabel("Temperature (K)")
    for ax in [new_ax, old_ax]:
        ax.set_xlabel(r'Radius $R_{\bigoplus}$')
        ax.grid(alpha=0.4)
    fig.save("iron_layer_particles.png", format='png')


def thermal_profile(new_particles, old_particles, earth_radius=6371 * 1000):
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.0, "wspace": 0.10})
    fig.patch.set_facecolor('xkcd:black')
    new_ax, old_ax = axs.flatten()[0], axs.flatten()[1]
    new_ax.set_title("New EoS")
    old_ax.set_title("Old EoS")
    new_ax.set_ylabel("Temperature (K)")
    for ax in [new_ax, old_ax]:
        ax.set_xlabel(r'Radius $R_{\bigoplus}$')
        ax.grid(alpha=0.4)
    for label in labels.keys():
        new_ax.scatter(
            [p.distance / earth_radius for p in new_particles if p.label == "PLANET" and p.tag == label and p.distance / (6371 * 1000) < 1.5],
            [p.temperature for p in new_particles if p.label == "PLANET" and p.tag == label and p.distance / (6371 * 1000) < 1.5],
            label=labels[label],
            s=1
        )
        old_ax.scatter(
            [p.distance / earth_radius for p in old_particles if p.label == "PLANET" and p.tag == label and p.distance / (6371 * 1000) < 1.5],
            [p.temperature for p in old_particles if p.label == "PLANET" and p.tag == label and p.distance / (6371 * 1000) < 1.5],
            label=labels[label],
            s=1
        )
    legend = new_ax.legend()
    for handle in legend.legendHandles:
        try:
            handle.set_sizes([3.0])
        except:
            pass

    plt.savefig("thermal_profile.png", format='png')


new_particles, new_time = get_particles(path=new_path, number_processes=new_processes, time=time,
                                        solve=False, formatted=True)
old_particles, old_time = get_particles(path=old_path, number_processes=old_processes, time=time,
                                        solve=False, formatted=True)
new_particles_t0, new_time_t0 = get_particles(path=new_path_unformatted, number_processes=new_processes, time=start_time,
                                        solve=True, formatted=False)
old_particles_t0, old_time_t0 = get_particles(path=old_path_unformatted, number_processes=old_processes, time=start_time,
                                        solve=True, formatted=False)

new_t0_dict = map_t0_particles(t_0_particles=new_particles_t0)
old_t0_dict = map_t0_particles(t_0_particles=old_particles_t0)
new_delta_t = iron_layer_delta_temperature(t_0_particles=new_particles_t0, t_particles=new_particles)
old_delta_t = iron_layer_delta_temperature(t_0_particles=old_particles_t0, t_particles=old_particles)
fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.0, "wspace": 0.10})
fig.patch.set_facecolor('xkcd:black')
new_ax, old_ax = axs.flatten()[0], axs.flatten()[1]
new_ax.scatter(
    [p.distance / (6371 * 1000) for p in new_particles if p.distance / (6371 * 1000) < 1.5],
    new_delta_t,
    color='magenta',
    s=2
)
old_ax.scatter(
    [p.distance / (6371 * 1000) for p in old_particles if p.distance / (6371 * 1000) < 1.5],
    old_delta_t,
    color='magenta',
    s=2
)
new_ax.set_title("New EoS")
old_ax.set_title("Old EoS")
new_ax.set_ylabel("Delta T (K)")
for ax in [new_ax, old_ax]:
    ax.set_xlabel(r'Radius $R_{\bigoplus}$')
    ax.grid(alpha=0.4)
plt.savefig("delta_T_iron_layer.png", format='png')


new_delta_rho = iron_layer_density_delta(t_0_particles=new_particles_t0, t_particles=new_particles)
old_delta_rho = iron_layer_density_delta(t_0_particles=old_particles_t0, t_particles=old_particles)
fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.0, "wspace": 0.10})
fig.patch.set_facecolor('xkcd:black')
new_ax, old_ax = axs.flatten()[0], axs.flatten()[1]
new_ax.scatter(
    [p.distance / (6371 * 1000) for p in new_particles if p.distance / (6371 * 1000) < 1.5 and p.tag == 3],
    new_delta_rho,
    color='magenta',
    s=2
)
old_ax.scatter(
    [p.distance / (6371 * 1000) for p in old_particles if p.distance / (6371 * 1000) < 1.5 and p.tag == 3],
    old_delta_rho,
    color='magenta',
    s=2
)
new_ax.set_title("New EoS")
old_ax.set_title("Old EoS")
new_ax.set_ylabel("Delta Density (kg/m3)")
for ax in [new_ax, old_ax]:
    ax.set_xlabel(r'Radius $R_{\bigoplus}$')
    ax.grid(alpha=0.4)
plt.savefig("delta_rho_iron_layer.png", format='png')




thermal_profile(new_particles=new_particles, old_particles=old_particles)
