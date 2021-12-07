#!/usr/bin/env python3
import pandas as pd
import numpy as np
from math import sqrt
from statistics import mean
import matplotlib.pyplot as plt

new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
start_time = 0
end_time = 2500
increment = 100

def calculate_smth_lengths(particle_density, particle_mass, dimension=3, eta=1.2):
    """
    Returns a list of smoothing lengths h.
    :param particles:
    :param dimension:
    :param eta:
    :return:
    """
    return [
        eta * (particle_mass / particle_density) ** (1 / dimension)
    ]

def calculate_w_ij(particle_i_position, particle_j_position, particle_i_velocity, particle_j_velocity):
    dr = [i - j for i, j in zip(particle_i_position, particle_j_position)]
    dv = [i - j for i, j in zip(particle_i_velocity, particle_j_velocity)]
    norm_dr = np.linalg.norm(dr)
    norm_dv = np.linalg.norm(dv)
    if norm_dr * norm_dv < 0:
        return (norm_dr * norm_dv) / sqrt(norm_dr ** 2)
    else:
        return 0

def calculate_v_sig(w_ij, particle_i_soundspeed, particle_j_soundspeed):
    return (particle_i_soundspeed + particle_j_soundspeed) - (3 * w_ij)

def calculate_vsig_max(curr_vsig, candidate_vsig):
    if curr_vsig < candidate_vsig:
        return candidate_vsig
    return curr_vsig


