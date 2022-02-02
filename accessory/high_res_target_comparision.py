import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("dark_background")

low_res_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/5_new/tar.dat"
high_res_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/5_new_high/target_5_new_high"
y_param = "density"

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)

low_df = pd.read_csv(low_res_path, skiprows=2, header=None, delimiter="\t")
low_df['radius'] = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) / (6371 * 1000) for i, j, k in zip(low_df[3], low_df[4], low_df[5])]
low_silicate = low_df[low_df[1] == 0]
low_iron = low_df[low_df[1] == 1]

high_df = pd.read_csv(high_res_path, skiprows=2, header=None, delimiter="\t")
high_df['radius'] = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) / (6371 * 1000) for i, j, k in zip(high_df[3], high_df[4], high_df[5])]
high_silicate = high_df[high_df[1] == 0]
high_iron = high_df[high_df[1] == 1]

ax.scatter(
    low_silicate['radius'], low_silicate[y_param], color='blue', s=2, label="Low Res Silicate)"
)
ax.scatter(
    low_iron['radius'], low_iron[y_param], color='blue', s=2, label="Low Res Iron)"
)
ax.scatter(
    high_silicate['radius'], high_silicate[y_param], color='magenta', s=2, label="High Res Silicate)"
)
ax.scatter(
    high_iron['radius'], high_iron[y_param], color='green', s=2, label="High Res Iron)"
)

ax.set_xlabel(r"Radius ($R_{\bigoplus}$)")
ax.set_ylabel(y_param.capitalize())
ax.grid(alpha=0.4)
ax.legend(loc='upper left')
ax.title("Low vs. High Res Target (5 kg/m3 Cutoff Denisty)")
plt.savefig("low_vs_high_comparision.png", format='png', dpi=200)
