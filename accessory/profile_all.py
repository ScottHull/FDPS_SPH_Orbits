#!/usr/bin/env python3
import os
import shutil
from random import randint
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

start_time = 80
end_time = 3000
interval = 20
number_processes = 200
min_norm = 0
max_norm = 10000
path = "/home/theia/scotthull/gi_new_eos"
eos = "src/phase_data/forst_STS.rho_u.txt"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_contour_GI"