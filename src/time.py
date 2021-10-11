import os
import csv
from src.combine import CombineFile

def get_all_iterations_and_times(min_iteration, max_iteration, number_processes, path):
    sampled_times = []
    fname_template = path + "/results.{}_{}_{}.dat"
    for time in range(min_iteration, max_iteration + 1):
        fname = fname_template.format(str(time).zfill(5), str(number_processes).zfill(5),
                                      str(number_processes - 1).zfill(5))
        with open(fname, 'r') as infile:
            reader = csv.reader(infile, delimiter="\t")
            sim_time = float(list(next(reader))[0])
            sampled_times.append((time, sim_time))
            infile.close()
    return sampled_times

def get_nearest_iteration_to_time(time, sampled_times):
    differences = []
    for iteration, time_candidate in sampled_times:
        diff = abs(time - time_candidate)
        differences.append((iteration, diff))
    min_diff = min([i[1] for i in differences])
    min_diff_index = [i[1] for i in differences].index(min_diff)
    return differences[min_diff_index][0]

def seconds_to_hours(seconds):
    # 60 seconds in 1 minute, 60 minutes in 1 hour
    return seconds * (1.0 / 60.0) * (1.0 / 60.0)

def match_particle_properties_between_iterations(particles1, particles2, property):
    particles = []
    for p1 in particles1:
        for p2  in particles2:
            if p1.particle_id == p2.particle_id:
                particles.append([p1, p2])
                break
    return particles
