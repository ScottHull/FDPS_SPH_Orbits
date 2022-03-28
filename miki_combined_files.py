import os
import csv
import pandas as pd
import numpy as np
from random import randint
import multiprocessing as mp

from src.combine import CombineFile

path_to_unformatted = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new_high/5_b073_new_high"
to_path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new_high/miki_combined_format"
min_iteration = 0
max_iteration = 1300
increment = 5
number_processes = 500

header = ["X", "Y", "Z", "VX", "VY", "VZ", "RADIUS", "ID", "TYPE", "VALUE1"]

if not os.path.exists(to_path):
    os.mkdir(to_path)

def make_file(iteration):
    to_fname = "merged_{}_{}.dat".format(max_iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=max_iteration,
                     output_path=path_to_unformatted, to_fname=to_fname)
    combined_file = cf.combine()
    formatted_time = round(cf.sim_time * 0.000277778, 2)
    particle_count = cf.total_particles
    f = os.getcwd() + "/{}".format(to_fname)
    df = pd.read_csv(f, skiprows=2, header=None, delimiter='\t')
    os.remove(f)
    df_new = pd.DataFrame({})
    df_new['X'], df_new['Y'], df_new['Z'] = df[3], df[4], df[5]
    df_new['VX'], df_new['VY'], df_new['VZ'] = df[6], df[7], df[8]
    df_new['RADIUS'] = [((i ** 2 + j ** 2 + k ** 2) ** (1 / 2)) / (10 ** 6) for i, j, k in zip(df[3], df[4], df[5])]
    df_new['ID'] = df[0]
    df_new['TYPE'] = [0 if i % 2 == 0 else 1 for i in df[2]]
    df_new['VALUE1'] = df[13]
    new_fname = to_path + "/{}.dat".format(iteration)
    df_new.to_csv(new_fname, sep="\t")
    with open(new_fname, 'r+') as infile:
        content = infile.read()
        infile.seek(0, 0)
        infile.write(formatted_time + "\n" + particle_count + "\n" + content)
    infile.close()

pool = mp.Pool(5)
pool.map(make_file, [[iteration] for iteration in
                     np.arange(min_iteration, max_iteration + increment, increment)])
pool.close()
pool.join()
