import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

fname = "/home/theia/scotthull/examine_fdps_sph_outputs/merged_142.dat"
radius_planet = 3800 * 1000

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature", "smoothing length"]
df = pd.read_csv(fname, skiprows=2, header=None, delimiter="\t", names=headers, index_col='id')
df['radius'] = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in zip(df['x'], df['y'], df['z'])]
df_silicate = df[df['tag'] % 2 == 0]

interp_functions = {
    header: interp1d(df_silicate['radius'], df_silicate[header]) for header in headers[9:]
}

# get all particles outside the planet
df_ejected = df_silicate[df_silicate['radius'] > radius_planet]
num_particles_to_fix = len(df_ejected)
print(f"there are {num_particles_to_fix} particles outside the planet)")

count = 0
# place these particles in random locations within the planet
for particle in df_ejected.index.tolist():
    # delete the particle
    df = df.drop(particle)

    # count += 1
    # print(f"particle {count} of {num_particles_to_fix} (particle id: {particle})")
    # # get random radius within planet
    # radius = np.random.uniform(min(df_silicate['radius']), radius_planet - (800 * 1000))
    # # get random theta
    # theta = np.random.uniform(0, 2 * np.pi)
    # # get random phi
    # phi = np.random.uniform(0, np.pi)
    # # convert to cartesian
    # x = radius * np.sin(phi) * np.cos(theta)
    # y = radius * np.sin(phi) * np.sin(theta)
    # z = radius * np.cos(phi)
    # print(f"old position: {df.at[particle, 'x']}, {df.at[particle, 'y']}, {df.at[particle, 'z']}...radius: {df.at[particle, 'radius'] / 1000} km")
    # print(f"new position: {x}, {y}, {z}...radius: {radius / 1000} km")
    # # update dataframe
    # df.at[particle, 'x'] = x
    # df.at[particle, 'y'] = y
    # df.at[particle, 'z'] = z
    # df.at[particle, 'radius'] = radius
    # df.at[particle, 'vx'] = 0
    # df.at[particle, 'vy'] = 0
    # df.at[particle, 'vz'] = 0
    # # interpolate the rest of the headrs as a function of radius
    # for header in headers[9:-1]:
    #     # print(f'interpolating {header}')
    #     f = interp_functions[header]
    #     df.at[particle, header] = f(radius)




# # get all particles in df whose ids are in df_ejected
# df_fixed = df[df.index.isin(df_ejected.index.tolist())]
# # plot to check
# fig, axs = plt.subplots(2, 2, figsize=(10, 5))
# axs = axs.flatten()
# for index, header in enumerate(['density', 'internal energy', 'entropy', 'temperature']):
#     # scatter the particles and color by 'tag' column
#     axs[index].scatter(
#         df[df['tag'] % 2 == 0]['radius'] / 1000, df[df['tag'] % 2 == 0][header], s=2, alpha=0.5, color='red', label="Silicate Particles"
#     )
#     axs[index].scatter(
#         df[df['tag'] % 2 == 1]['radius'] / 1000, df[df['tag'] % 2 == 1][header], s=2, alpha=0.5, color='blue', label="Iron Particles"
#     )
#     axs[index].scatter(
#         df_fixed['radius'] / 1000, df_fixed[header], s=2, alpha=0.5, color='green', label="Fixed Particles"
#     )
#     # axs[index].scatter(
#     #     df['radius'] / 1000, df[header], s=2, alpha=0.5, label="All Particles"
#     # )
#     # plot the radius of the planet
#     axs[index].axvline(radius_planet / 1000, color='k', linestyle='--', label="Planet Radius")
#     # axs[index].set_title(header)
#     axs[index].set_ylabel(header)
# axs[0].legend(loc='upper right')
# for ax in axs:
#     ax.grid()




# remove the smoothing length column
df = df.drop(columns=['smoothing length', 'radius'])
# check for nan values in the dataframe and print out their ids if they exist
nan_values = df[df.isna().any(axis=1)]
if len(nan_values) > 0:
    print("there are nan values in the dataframe")
    print(nan_values.index.tolist())
    # print the corresponding rows
    print(df.loc[nan_values.index.tolist()])
    # delete these rows
    df = df.drop(nan_values.index.tolist())
# save the dataframe
df.to_csv('body.dat', sep='\t', index=True, header=False)
# plt.savefig('test.png', format='png')
