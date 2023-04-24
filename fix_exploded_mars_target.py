import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

fname = "merged_142.dat"
radius_planet = 1000

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature", "smoothing length"]
df = pd.read_csv(fname, skiprows=2, header=None, delimiter="\t", names=headers)
df['radius'] = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in zip(df['x'], df['y'], df['z'])]
df_silicate = df[df['tag'] % 2 == 0]

# get all particles outside the planet
df_ejected = df_silicate[df_silicate['x'] > radius_planet]

# place these particles in random locations within the planet
for particle in df_ejected.index.tolist():
    # get random radius within planet
    radius = np.random.uniform(0, radius_planet)
    # get random theta
    theta = np.random.uniform(0, 2 * np.pi)
    # get random phi
    phi = np.random.uniform(0, np.pi)
    # convert to cartesian
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)
    # update dataframe
    df.at[particle, 'x'] = x
    df.at[particle, 'y'] = y
    df.at[particle, 'z'] = z
    df.at[particle, 'vx'] = 0
    df.at[particle, 'vy'] = 0
    df.at[particle, 'vz'] = 0
    # interpolate the rest of the headrs as a function of radius
    for header in headers[9:]:
        f = interp1d(df_silicate['radius'], df_silicate[header])
        df.at[particle, header] = f(radius)

# plot to check
fig, axs = plt.subplots(2, 2, figsize=(10, 5))
axs = axs.flatten()
for index, header in enumerate(['density', 'internal energy', 'entropy', 'temperature']):
    axs[index].scatter(
        df['radius'] / 1000, df[header], s=2, alpha=0.5, label="All Particles"
    )
    # axs[index].set_title(header)
    axs[index].set_ylabel(header)
for ax in axs:
    ax.grid()

plt.savefig('test.png', format='png')
