from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

infile = "/Users/scotthull/Desktop/tar_old_eos.dat"
df = pd.read_csv(infile, skiprows=2, header=None, delimiter="\t")

particle_id, tag, position, density, internal_energy, pressure, entropy, temperature = df[0], df[1], zip(df[3], df[4],
                                                                                                         df[5]), df[9], \
                                                                                       df[10], df[11], df[13], df[14]
radius = [sqrt(x ** 2 + y ** 2 + z ** 2) for x, y, z in position]
df['radius'] = radius
df['tag'] = tag
to_plot = [("internal_energy", internal_energy), ("density", density), ("pressure", pressure), ("entropy", entropy),
           ("temperature", temperature)]
for i, j in to_plot:
    df[i] = j

sns.set()
fig, ax = plt.subplots(3, 2, figsize=(8, 8))
for index, data in enumerate(to_plot):
    label = data[0]
    sns.scatterplot(data=df, x=df['radius'], y=label, ax=ax.flatten()[index], hue="tag")

plt.suptitle("Target (Old EoS)")
plt.tight_layout()
plt.show()
