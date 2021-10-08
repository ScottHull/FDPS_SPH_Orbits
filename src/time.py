import os
from src.combine import CombineFile

def get_nearest_iteration_to_time(time, min_iteration, max_iteration, number_processes, path):
    sampled_times = []
    for i in range(min_iteration, max_iteration + 1):
        to_fname = "tmp.dat"
        cf = CombineFile(num_processes=number_processes, time=time, output_path=path, to_fname=to_fname)
        t = cf.sim_time
        sampled_times.append(t)
        if i == max_iteration:
            return i
        if abs(sampled_times[-2] - t) < abs(sampled_times[-1] - t):
            return i - 1

def seconds_to_hours(seconds):
    # 60 seconds in 1 minute, 60 minutes in 1 hour
    return seconds * (1.0 / 60.0) * (1.0 / 60.0)
