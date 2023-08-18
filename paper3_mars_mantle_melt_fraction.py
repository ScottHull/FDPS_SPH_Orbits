#!/usr/bin/env python3
import os
import shutil
import csv
import pandas as pd
import numpy as np
import string
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.animate import animate
from src.combine import CombineFile


plt.rcParams.update({'font.size': 14, })
plt.style.use("dark_background")
# plt.style.use('seaborn-pastel')

base_path = "/home/theia/scotthull/Paper2_SPH/gi/"
num_processes = 600
paths = [['500_mars', "Mars " + r"($b=0.73$)"]]
initial_iteration = 0
post_impact_iteration = 20

initial_file = CombineFile(num_processes=num_processes, time=initial_iteration, output_path=f"{base_path}{paths[0][0]}/{paths[0][0]}")
initial_df = initial_file.combine_df()
post_impact_file = CombineFile(num_processes=num_processes, time=post_impact_iteration, output_path=f"{base_path}{paths[0][0]}/{paths[0][0]}")
post_impact_df = post_impact_file.combine_df()

delta_S = {i: post_impact_df.loc[i, 'id']['entropy'] - initial_df.loc[i, 'id']['entropy'] for i in initial_df['index'].tolist()}
print(len([i for i in delta_S.values() if i >= 500]) / len(delta_S.values())
