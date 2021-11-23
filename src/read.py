import csv
import pandas as pd

from src.identify import ParticleMap

def get_particles_from_formatted(path, iteration, solve=False):
    f = "{}/{}.csv".format(path, iteration)
    time = 0
    num_particles = 0
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter=",")
        time = float(next(reader)[0])
        num_particles = float(next(reader)[0])
        infile.close()
    df = pd.read_csv(f, delimiter=",", skiprows=2, index_col="id")
    return time, num_particles, df
