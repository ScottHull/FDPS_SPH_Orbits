#!/usr/bin/env python3
import os
from random import randint

from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import write_report_at_time

iteration = 1500
number_processes = 200
path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new/5_b073_new"
output_name = path.split("/")[-1]

if "new" in path:
    phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
else:
    phase_path = "src/phase_data/duniteN_vapour_curve.txt"

to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
combined_file = cf.combine()
formatted_time = cf.sim_time
f = os.getcwd() + "/{}".format(to_fname)
pm = ParticleMap(path=f, center=True, relative_velocity=False)
particles = pm.collect_particles(find_orbital_elements=True)
os.remove(to_fname)
pm.solve(particles=particles, phase_path=phase_path, report_name="{}-report.txt".format(output_name),
         iteration=iteration, simulation_time=formatted_time)
write_report_at_time(particles=particles, fname="{}.csv".format(output_name))
