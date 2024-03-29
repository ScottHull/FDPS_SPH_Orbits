import os
from random import randint
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from src.combine import CombineFile

# use the dark background style
plt.style.use("dark_background")

path = "/home/theia/scotthull/Paper1_SPH/gi/500_b073_new/500_b073_new"
iteration = 30
number_processes = 200
square_scale = 2e7

to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
combined_file = cf.combine()
formatted_time = round(cf.sim_time * 0.000277778, 2)
f = os.getcwd() + "/{}".format(to_fname)
df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t",
                 names=["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                        "potential energy", "entropy", "temperature"])
# sort the df by z
df = df.sort_values(by=['z'])
os.remove(f)

# get a reds color map
normalizer = Normalize(2000, 10000)
cmap = cm.get_cmap('Reds').reversed()
# scatter the points with x and y as the axes, and the color mapped to temperature
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.scatter(df['x'], df['y'], c=df['temperature'], cmap=cmap, norm=normalizer, alpha=1, s=0.3)
ax.set_xlim(-square_scale, square_scale)
# get rid of the axes
ax.set_axis_off()
# get rid of the border
ax.set_frame_on(False)

# save the figure
plt.tight_layout()
fig.savefig("test.png", dpi=300, bbox_inches='tight', pad_inches=0)
