#!/usr/bin/env python3
import os
from random import randint


from src.identify import ParticleMap
from src.combine import CombineFile
from src.vapor import calc_vapor_mass_fraction

end_iteration = 3000
number_processes = 200
output_name = "gi_new_eos_b_073"
path = '/home/theia/scotthull/1M/gi_new_eos_b_073'

if "new" in path:
    phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
else:
    phase_path = "src/phase_data/duniteN_vapour_curve.txt"

if "formatted" in path:
    f = path + "/{}.csv".format(end_iteration)
    pm = ParticleMap(path=f, center=True, relative_velocity=False, formatted=True)
    particles = pm.collect_particles(find_orbital_elements=False)
    vmf = calc_vapor_mass_fraction(particles=particles, phase_path=phase_path, only_disk=True) * 100.0

    pm.profile_report.update({"disk vmf": vmf})
    pm.profile_report.update({"name": output_name})
    pm.profile_report.update({'phase_path': phase_path})
    pm.profile_report.update({'iteration profiled (corresponds to output number)': end_iteration})
    pm.report(name=output_name)
else:
    to_fname = "merged_{}_{}.dat".format(end_iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=end_iteration, output_path=path, to_fname=to_fname)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    f = os.getcwd() + "/{}".format(to_fname)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=True)
    pm.solve(particles=particles, phase_path=phase_path, report_name="{}.txt".format(output_name), iteration=end_iteration, simulation_time=formatted_time)
