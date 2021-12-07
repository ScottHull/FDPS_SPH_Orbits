#!/usr/bin/env python3
import pandas as pd
import numpy as np
from math import sqrt
from statistics import mean
import matplotlib.pyplot as plt

from src.interpolation import GenericTrilinearInterpolation

time = 1000
mean_artificial_viscosity = 1.5
new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
new_eos_silicate_table = "src/phase_data/forst_STS.rho_u.txt"
new_eos_iron_table = "src/phase_data/iron_STS.rho_u.txt"
old_eos_silicate_table = "src/phase_data/duniteN.rho_u.txt"
old_eos_iron_table = "src/phase_data/ironC.rho_u.txt"

new_silicate_eos = pd.read_fwf(new_eos_silicate_table, skiprows=2,
                               names=["density", "internal_energy", "temperature", "pressure",
                                      "soundspeed", "entropy"])
new_iron_eos = pd.read_fwf(new_eos_iron_table, skiprows=2,
                           names=["density", "internal_energy", "temperature", "pressure",
                                  "soundspeed", "entropy"])
old_silicate_eos = pd.read_fwf(old_eos_silicate_table, skiprows=2,
                               names=["density", "internal_energy", "temperature", "pressure",
                                      "soundspeed", "entropy"])
old_iron_eos = pd.read_fwf(old_eos_iron_table, skiprows=2,
                           names=["density", "internal_energy", "temperature", "pressure",
                                  "soundspeed", "entropy"])

new_interp_silicate = GenericTrilinearInterpolation(
    var1_array=list(new_silicate_eos['density']), var2_array=list(new_silicate_eos['internal_energy']),
    var3_array=list(new_silicate_eos['soundspeed'])
)
new_interp_iron = GenericTrilinearInterpolation(
    var1_array=list(new_iron_eos['density']), var2_array=list(new_iron_eos['internal_energy']),
    var3_array=list(new_iron_eos['soundspeed'])
)
old_interp_silicate = GenericTrilinearInterpolation(
    var1_array=list(old_silicate_eos['density']), var2_array=list(old_silicate_eos['internal_energy']),
    var3_array=list(old_silicate_eos['soundspeed'])
)
old_interp_iron = GenericTrilinearInterpolation(
    var1_array=list(old_iron_eos['density']), var2_array=list(old_iron_eos['internal_energy']),
    var3_array=list(old_iron_eos['soundspeed'])
)


def calculate_smoothing_length(mass_particle, density_particle):
    return 1.2 * ((mass_particle / density_particle) ** (1 / 3))


def velocity_sph(mean_av, soundspeed, mass_particle, density_particle):
    smoothing_length = calculate_smoothing_length(mass_particle, density_particle)
    return (soundspeed * smoothing_length * mean_av) / 8.0


def numerical_viscosity_time(particle_radius, mean_av, mass_particle, density_particle, internal_energy_particle,
                             interp_class):
    soundspeed = interp_class.interpolate(var1=density_particle, var2=internal_energy_particle)
    v_sph = velocity_sph(mean_av, soundspeed, mass_particle, density_particle)
    return ((particle_radius ** 2) / v_sph) * 0.000277778  # seconds --> hours


print("At time : {}".format(time))
new_f, old_f = new_path + "/{}.csv".format(time), old_path + "/{}.csv".format(time)
new_file = pd.read_csv(new_f, skiprows=2).to_dict('list')
old_file = pd.read_csv(old_f, skiprows=2).to_dict('list')

new_tau_silicate = [
    [   
        radius,
        numerical_viscosity_time(
            radius,
            mean_artificial_viscosity,
            new_file['mass'][index],
            new_file['density'][index],
            new_file['internal_energy'][index],
            new_silicate_eos
        )
    ]
    for index, radius in enumerate(new_file['radius']) if new_file['tag'][index] % 2 == 0 
                                                          and new_file['label'][index] == "DISK"
]
old_tau_silicate = [
    [   
        radius,
        numerical_viscosity_time(
            radius,
            mean_artificial_viscosity,
            old_file['mass'][index],
            old_file['density'][index],
            old_file['internal_energy'][index],
            old_silicate_eos
        )
    ]
    for index, radius in enumerate(old_file['radius']) if old_file['tag'][index] % 2 == 0 
                                                          and old_file['label'][index] == "DISK"
]
# new_tau_iron = [
#     [
#         radius,
#         numerical_viscosity_time(
#             radius,
#             mean_artificial_viscosity,
#             new_file['mass'][index],
#             new_file['density'][index],
#             new_file['internal_energy'][index],
#             new_iron_eos
#         )
#     ]
#     for index, radius in enumerate(new_file['radius']) if new_file['tag'][index] % 2 == 0
#                                                           and new_file['label'][index] == "DISK"
# ]
# old_tau_iron = [
#     [
#         radius,
#         numerical_viscosity_time(
#             radius,
#             mean_artificial_viscosity,
#             old_file['mass'][index],
#             old_file['density'][index],
#             old_file['internal_energy'][index],
#             old_iron_eos
#         )
#     ]
#     for index, radius in enumerate(old_file['radius']) if old_file['tag'][index] % 2 == 0
#                                                           and old_file['label'][index] == "DISK"
# ]


plt.style.use("dark_background")
fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.14})
fig.patch.set_facecolor('xkcd:black')

new_fig, old_fig = axs.flatten()[0], axs.flatten()[1]

new_fig.scatter(
    [i[0] / (6371 * 1000) for i in new_tau_silicate],
    [i[1] for i in new_tau_silicate],
    s=2
)
old_fig.scatter(
    [i[0] / (6371 * 1000) for i in old_tau_silicate],
    [i[1] for i in old_tau_silicate],
    s=2
)
for ax in axs.flatten():
    ax.set_xlabel(r'Radius $R_{\bigoplus}$')
    ax.grid(alpha=0.4)
new_fig.set_title("Numerical Viscosity (New EoS)")
old_fig.set_title("Numerical Viscosity (Old EoS)")
new_fig.set_ylabel("tau_sph")

plt.savefig("numerical_viscosity.png", format='png')
