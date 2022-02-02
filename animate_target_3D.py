import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.combine import CombineFile
from src.animate import animate

plt.style.use("dark_background")

from_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/5_new_high/target_5_new_high"
to_path = "animate_target"
start_iteration = 0
end_iteration = 60
num_proc = 400
square_scale = 5.5

def get_time(f, local=True):
    formatted_time = None
    if local:  # if not reading from remote server
        with open(f, 'r') as infile:
            reader = csv.reader(infile, delimiter="\t")
            formatted_time = float(next(reader)[0])
        infile.close()
    else:
        formatted_time = float(next(f))
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours

if not os.path.exists(to_path):
    os.mkdir(to_path)

for iteration in np.arange(start_iteration, end_iteration + 1, 1):
    to_fname = "target_{}.csv".format(iteration)
    c = CombineFile(num_processes=num_proc, time=iteration, output_path=from_path, to_fname=to_fname).combine()
    df = pd.read_csv(to_fname, skiprows=2, header=None, delimiter="\t")
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        df[3], df[4], df[5], s=2
    )
    ax.set_xlabel("x")
    ax.set_xlabel("y")
    ax.set_xlabel("z")
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_zlim(-square_scale, square_scale)
    ax.set_title("Target ({} hrs)".format(get_time(to_fname)))
    plt.savefig(to_path + "/{}.png".format(iteration), format='png', dpi=200)
    os.remove(to_fname)

animate(
    start_time=start_iteration,
    end_time=end_iteration,
    interval=1,
    path=to_path,
    fps=5,
    filename="target_animate.mp4",
)
