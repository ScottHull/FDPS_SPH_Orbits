#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.new_and_old_eos import seconds_to_hours

start_time = 0
end_time = 3000
interval = 20
number_processes = 200
path = "/home/theia/scotthull/gi_new_eos"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_contour_GI"

for o in [output]:
    if os.path.exists(o):
        shutil.rmtree(o)
    os.mkdir(o)




