#!/usr/bin/env python3
import pandas as pd
import numpy as np
from statistics import mean
import matplotlib.pyplot as plt

new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
start_time = 0
end_time = 2500
increment = 100

def calculate_smth_lengths(particles, dimension=3, eta=1.2):
    """
    Returns a list of smoothing lengths h.
    :param particles:
    :param dimension:
    :param eta:
    :return:
    """
    return [
        eta * (mass / particles["density"][index]) ** (1 / dimension) for mass, index in particles['mass']
    ]

