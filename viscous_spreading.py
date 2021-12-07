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
new_eos_silicate_table = "src/phase_data/forstSTS.rho_u.txt"
new_eos_iron_table = "src/phase_data/ironSTS.rho_u.txt"
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
    var1_array=new_silicate_eos['density'], var2_array=new_silicate_eos['internal_energy'],
    var3_array=new_silicate_eos['soundspeed']
)
new_interp_iron = GenericTrilinearInterpolation(
    var1_array=new_iron_eos['density'], var2_array=new_iron_eos['internal_energy'],
    var3_array=new_iron_eos['soundspeed']
)
old_interp_silicate = GenericTrilinearInterpolation(
    var1_array=old_silicate_eos['density'], var2_array=old_silicate_eos['internal_energy'],
    var3_array=old_silicate_eos['soundspeed']
)
old_interp_iron = GenericTrilinearInterpolation(
    var1_array=old_iron_eos['density'], var2_array=old_iron_eos['internal_energy'],
    var3_array=old_iron_eos['soundspeed']
)

print("At time : {}".format(time))
new_f, old_f = new_path + "/{}.csv".format(time), old_path + "/{}.csv".format(time)
new_file = pd.read_csv(new_f, skiprows=2).to_dict('list')
old_file = pd.read_csv(old_f, skiprows=2).to_dict('list')

new_soundspeeds_silicate = [new_interp_silicate.interpolate(var1=density, var2=new_file['internal_energy'][index]) for
                            density, index in enumerate(new_file['density']) if
                            new_file['label'][index] == "DISK" and new_file['tag'][index] % 2 == 0]
new_soundspeeds_iron = [new_interp_iron.interpolate(var1=density, var2=new_file['internal_energy'][index]) for
                            density, index in enumerate(new_file['density']) if
                            new_file['label'][index] == "DISK" and new_file['tag'][index] % 2 == 0]
old_soundspeeds_silicate = [old_interp_silicate.interpolate(var1=density, var2=old_file['internal_energy'][index]) for
                            density, index in enumerate(old_file['density']) if
                            old_file['label'][index] == "DISK" and old_file['tag'][index] % 2 == 0]
old_soundspeeds_iron = [old_interp_iron.interpolate(var1=density, var2=old_file['internal_energy'][index]) for
                            density, index in enumerate(old_file['density']) if
                            old_file['label'][index] == "DISK" and old_file['tag'][index] % 2 == 0]


