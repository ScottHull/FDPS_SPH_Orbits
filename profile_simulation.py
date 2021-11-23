import os
from random import randint


from src.identify import ParticleMap
from src.combine import CombineFile

end_iteration = 3000
number_processes = 200
output_name = "gi_new_eos_b_073"
path = '/home/theia/scotthull/1M/gi_new_eos_b_073'

if "new" in path:
    phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
else:
    phase_path = "src/phase_data/duniteN_vapour_curve.txt"

to_fname = "merged_{}_{}.dat".format(end_iteration, randint(0, 100000))
cf = CombineFile(num_processes=number_processes, time=end_iteration, output_path=path, to_fname=to_fname)
combined_file = cf.combine()
formatted_time = cf.sim_time
f = os.getcwd() + "/{}".format(to_fname)
pm = ParticleMap(path=f, center=True, relative_velocity=False)
particles = pm.collect_particles(find_orbital_elements=False)
pm.solve(particles=particles, phase_path=phase_path, report_name="{}.txt".format(output_name))
