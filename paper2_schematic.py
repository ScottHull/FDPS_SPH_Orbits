import os
from random import randint
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from src.combine import CombineFile

path = "/home/theia/scothull/Paper1_SPH/gi/500_b073_new/500_b073_new/"
iteration = 20
number_processes = 200
square_scale = 2e7
file_format = "results.{}_{}_{}.dat"

to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
combined_file = cf.combine()
formatted_time = round(cf.sim_time * 0.000277778, 2)
f = os.getcwd() + "/{}".format(to_fname)
df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t",
                 names=["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                        "potential energy", "entropy", "temperature"])
os.remove(f)

# get a reds color map
normalizer = Normalize(2000, 10000)
cmap = cm.get_cmap('reds')
# scatter the points with x and y as the axes, and the color mapped to temperature
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, aspect='equal')
ax.scatter(df['x'], df['y'], c=df['temperature'], cmap=cmap, norm=normalizer, alpha=1, s=0.3)
ax.set_xlim(-square_scale, square_scale)
# get rid of the axes
ax.set_axis_off()
# get rid of the border
ax.set_frame_on(False)

# save the figure
fig.savefig("test.png", dpi=300, bbox_inches='tight', pad_inches=0)
