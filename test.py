import os
import shutil
import numpy as np
from math import cos
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

time = 3000
number_processes = 100
path = "/Users/scotthull/Desktop/merged_3000.dat"

pm = ParticleMap(path=path, center=True, relative_velocity=False)
particles = pm.collect_particles()
pm.solve(particles=particles)
