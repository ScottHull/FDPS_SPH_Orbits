import os
import pandas as pd
import numpy as np
from random import randint
import multiprocessing.dummy as mp

from src.combine import CombineFile
from src.identify import ParticleMap


class BuildReports:

    def __init__(self, to_dir, from_dir, start_time, end_time, number_processes, eos_phase_path, accessory_path, interval=1):
        self.to_dir = to_dir
        self.from_dir = from_dir
        self.start_time = start_time
        self.end_time = end_time
        self.number_processes = number_processes
        self.eos_phase_path = eos_phase_path
        self.interval = interval
        self.labels = {}
        self.accessory_path = accessory_path
        try:
            os.mkdir(accessory_path)
        except:
            pass
        self.__get_end_state()

    def __output_disk_state(self, particles, time, vmf):
        output = open("{}/{}.txt".format(self.accessory_path, time), 'w')
        avg_disk_entropy = [p.entropy for p in particles if p.label == "DISK"]
        avg_disk_entropy = sum(avg_disk_entropy) / len(avg_disk_entropy)
        disk_mass = sum([p.mass for p in particles if p.label == "DISK"])
        line = "time\t{}\ndisk vmf\t{}\navg disk entropy\t{}\ndisk mass\t{}\n".format(time, vmf, avg_disk_entropy, disk_mass)
        output.write(line)
        output.close()


    def __build_df_from_endstate(self, particles):
        try:
            end_labels = [self.labels[p.particle_id] for p in particles]
        except:
            end_labels = []
        return pd.DataFrame({
            "id": [p.particle_id for p in particles],
            "x": [p.position[0] for p in particles],
            "y": [p.position[1] for p in particles],
            "z": [p.position[2] for p in particles],
            "vx": [p.velocity[0] for p in particles],
            "vy": [p.velocity[1] for p in particles],
            "vz": [p.velocity[2] for p in particles],
            "mass": [p.mass for p in particles],
            "radius": [p.distance for p in particles],
            "tag": [p.tag for p in particles],
            "end_label": end_labels,
            "label": [p.label for p in particles],
            "entropy": [p.entropy for p in particles],
            "temperature": [p.temperature for p in particles],
            "density": [p.density for p in particles],
            "pressure": [p.pressure for p in particles],
            "internal_energy": [p.internal_energy for p in particles],
            "potential_energy": [p.potential_energy for p in particles],
            # "soundspeed": [p.soundspeed for p in particles]
        })

    def __get_end_state(self):
        particles = self.build_report_at_time(time=self.end_time, solve=True, save=False)
        self.labels = dict([(p.particle_id, p.label) for p in particles])
        self.build_report_at_time(time=self.end_time, solve=False, save=True)

    def build_report_at_time(self, time, pm=None, formatted_time=None, total_particles=None, solve=False, save=True):
        to_fname = "{}.csv".format(time)
        if not pm:
            fname = "merged_{}_{}_reportsproc.dat".format(time, randint(0, 1000000))
            cf = CombineFile(num_processes=self.number_processes, time=time, output_path=self.from_dir, to_fname=fname)
            combined_file = cf.combine()
            formatted_time = cf.sim_time
            total_particles = cf.total_particles
            f = os.getcwd() + "/{}".format(fname)
            pm = ParticleMap(path=f, center=True, relative_velocity=False)
            particles = pm.collect_particles(find_orbital_elements=solve)
            if solve:
                pm.solve(particles=particles, phase_path=self.eos_phase_path)
                self.__output_disk_state(time=time, particles=particles, vmf=pm.vmf)
            os.remove(f)
        else:
            particles = pm.particles
        if save:
            particle_df = self.__build_df_from_endstate(particles=particles)
            particle_df.to_csv(self.to_dir + "/{}".format(to_fname))
            with open(self.to_dir + "/{}".format(to_fname), 'r+') as infile:
                content = infile.read()
                infile.seek(0, 0)
                infile.write("{}\n{}\n".format(formatted_time, total_particles) + content)
            infile.close()
        return particles


    def make_reports(self, accessory_path, mp_pool_size=5):
        os.mkdir(accessory_path)
        pool = mp.Pool(mp_pool_size)
        for time in np.arange(self.start_time, self.end_time + self.interval, self.interval):
            pool.map(self.build_report_at_time, [time])
        pool.close()
        pool.join()
