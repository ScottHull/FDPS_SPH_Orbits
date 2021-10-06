from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

infile = "tar_new_eos.dat"
df = pd.read_csv(infile, skiprows=2, header=None)

particle_id, tag, position, density, internal_energy, pressure, entropy, temperature = df[0], df[1], zip(df[3], df[4],
                                                                                                         df[5]), df[9], \
                                                                                       df[10], df[11], df[13], df[14]
radius = [sqrt(x ** 2 + y ** 2 + z ** 2) for x, y, z in position]
df['radius'] = radius

sns.set()
g = sns.scatterplot(df=df, col="sex", hue="smoker")
g.map(sns.scatterplot, "total_bill", "tip", alpha=.7)
g.add_legend()
