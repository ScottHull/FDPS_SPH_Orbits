#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import string
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
# plt.style.use('seaborn-colorblind')

path = r"C:\Users\Scott\Desktop\b073_secondary_impact_struc_data.csv"
df = pd.read_csv(path, index_col="Unnamed: 0")
runs = df.keys()

fig, axs = plt.subplots(1, 2, figsize=(16, 9))
axs = axs.flatten()

fig.supxlabel("Runs")
for ax in axs:
    ax.grid(alpha=0.4)
axs[0].set_ylabel("Mass ($M_{L}$)")
axs[1].set_ylabel("Angular Momentum ($L_{EM}$)")

for i, label in zip(["MASS SI", "MASS TAIL"], ["Secondary Impactor", "Tail"]):
    axs[0].plot(
        range(0, len(runs)), [df[run][i] for run in runs], linewidth=2.0, label=label
    )
    axs[0].set_xticks(ticks=range(0, len(runs)))
    axs[0].set_xticklabels(runs, rotation=45, fontsize=10)

for i, label in zip(["ANGULAR MOMENTUM SI", "ANGULAR MOMENTUM TAIL"], ["Secondary Impactor", "Tail"]):
    axs[1].plot(
        range(0, len(runs)), [df[run][i] for run in runs], linewidth=2.0, label=label
    )
    axs[1].set_xticks(ticks=range(0, len(runs)))
    axs[1].set_xticklabels(runs, rotation=45, fontsize=10)

plt.xticks(rotation=45)

axs[0].legend()
plt.savefig("b073_secondary_impact_mass_and_am_mpl.png", format='png', dpi=200)
