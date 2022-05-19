from src import interpolation
import pandas as pd
import matplotlib.pyplot as plt

f = r"C:\Users\Scott\Desktop\imp_2000_new.dat"
headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
               "potential energy", "entropy", "temperature"]
imp_df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
imp_df = imp_df[imp_df['tag'] % 2 == 0]
imp_df['radius'] = ((imp_df['x'] ** 2 + imp_df['y'] ** 2 + imp_df['z'] ** 2) ** 0.5) / 1000
imp_df = imp_df[imp_df['radius'] < 3200]
df = pd.read_fwf("src/phase_data/forst_STS.rho_u.txt", sep='\t', header=None, skiprows=2)
rho_i, u_i, p_i, s_i = df[0], df[1], df[3], df[5]
i = interpolation.GenericTrilinearInterpolation(rho_i.tolist(), u_i.tolist(), p_i.tolist())
err_interpolated_p = [imp_df['pressure'][j] - i.interpolate(imp_df['density'][j], imp_df['internal energy'][j]) for j in imp_df.index]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    imp_df['radius'], err_interpolated_p, s=2, marker='.'
)

plt.show()
