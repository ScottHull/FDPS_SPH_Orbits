import os
import pandas as pd
import numpy as np
from random import randint
import multiprocessing.dummy as mp

from src.combine import CombineFile, make_particle_df
from src.identify import ParticleMap


class BuildReports:

    def __init__(self, to_dir, from_dir, start_time, end_time, number_processes, eos_phase_path, interval=1):
        self.to_dir = to_dir
        self.from_dir = from_dir
        self.start_time = start_time
        self.end_time = end_time
        self.number_processes = number_processes
        self.eos_phase_path = eos_phase_path
        self.interval = interval
        self.labels = {}

    def __build_df_from_endstate(self, particles):
        return pd.DataFrame({
            "id": [p.particle_id for p in particles],
            "position": [p.position for p in particles],
            "radius": [p.distance for p in particles],
            "tag": [p.tag for p in particles],
            "label": [self.labels[p.particle_id] for p in particles],
            "entropy": [p.entropy for p in particles],
            "temperature": [p.temperature for p in particles],
            "density": [p.density for p in particles],
            "pressure": [p.pressure for p in particles],
            "internal_energy": [p.internal_energy for p in particles],
        })

    def __get_end_state(self):
        particles = self.__build_report_at_time(time=self.end_time, solve=True)
        self.labels = dict([(p.particle_id, p.label) for p in particles])

    def __build_report_at_time(self, time, solve=False):
        fname = "merged_{}_{}.dat".format(time, randint(0, 1000000))
        cf = CombineFile(num_processes=self.number_processes, time=time, output_path=self.from_dir, to_fname=fname)
        combined_file = cf.combine()
        formatted_time = cf.sim_time
        total_particles = cf.total_particles
        f = os.getcwd() + "/{}".format(fname)
        pm = ParticleMap(path=f, center=True, relative_velocity=False)
        particles = pm.collect_particles(find_orbital_elements=False)
        if solve:
            pm.solve(particles=particles, phase_path=self.eos_phase_path)
        os.remove(f)
        particle_df = make_particle_df(particles=particles)
        particle_df.to_csv(self.to_dir + "/{}".format(fname))
        with open(self.to_dir + "/{}".format(fname), 'r+') as infile:
            content = infile.read()
            infile.seek(0, 0)
            infile.write("{}\n{}\n".format(formatted_time, total_particles) + content)
        infile.close()
        return particles

    def make_reports(self, mp_pool_size=5):
        self.__get_end_state()
        pool = mp.Pool(mp_pool_size)
        for time in np.arange(self.start_time, self.end_time + self.interval, self.interval):
            pool.map(self.__build_report_at_time, [time])
        pool.close()
        pool.join()
