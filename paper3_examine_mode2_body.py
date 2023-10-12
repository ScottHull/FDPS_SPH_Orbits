import numpy as np
import matplotlib.pyplot as plt

from src.combine import CombineFile

path = "/home/theia/scotthull/Paper3_SPH/tar-imp/target_mars_n_sph"
iteration = 71
num_processes = 600

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]

cf = CombineFile(
    num_processes=num_processes, time=iteration, output_path=path, to_fname=""
)
df = cf.combine_df()
df.columns = headers
df['radius'] = np.sqrt(df['x'] ** 2 + df['y'] ** 2 + df['z'] ** 2)

# first, make a 3D plot of the body in using xyz coordinates
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df['x'], df['y'], df['z'], marker=".", s=5)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_title("Mars (t={})".format(iteration))
plt.tight_layout()
plt.savefig(f"paper3_mode2_3D.png", dpi=200)

# then, make a 2x2 figure showing density, internal energy, pressure, and entropy as a function of radius
fig, axs = plt.subplots(2, 2, figsize=(10, 10))
axs = axs.flatten()
for ax, header in zip(axs, ["density", "internal energy", "pressure", "entropy"]):
    ax.scatter(df['radius'], df[header], marker=".", s=5)
    ax.set_xlabel("Radius (m)")
    ax.set_ylabel(header)
    ax.grid()
plt.tight_layout()
plt.savefig(f"paper3_mode2_scatter.png", dpi=200)
