from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt

from src.new_and_old_eos import get_particles

"""
We want to model the iron layer after the GI as a spherical shell.
Volume of a sphere: 4/3 pi r^3
Volume of a spherical shell: 4/3 pi (r_2 - r_1)^3
We need to calculate r_2 - r_1 (can't get it directly from simulation)
Know that density is rho = M / V, so the density of the shell is rho = M / [4/3 pi (r_2 - r_1)^3]
We can get M from the simulation, but what about density?
One way:
"""

time = 3000
new_path = "/home/theia/scotthull/1M/formatted_gi_new_eos_b_073"
old_path = "/home/theia/scotthull/1M/formatted_gi_old_eos_b_073"
new_processes = 200
old_processes = 200


plt.style.use("dark_background")

labels = {
    0: "Target Silicate",
    1: "Target Iron",
    2: "Impactor Silicate",
    3: "Impactor Iron"
}

def __mean_curve(x, y):
    combined = list(zip(x, y))
    sort = sorted(combined, key=lambda tup: tup[0])
    x = [i[0] for i in sort]
    y = [i[1] for i in sort]
    s = UnivariateSpline(x, y)
    xs = np.linspace(0, int(max(y)), 10000)
    ys = s(xs)
    return xs, ys

def thermal_profile(new_particles, old_particles, earth_radius = 6371 * 1000):

    fig, axs = plt.subplots(1, 2, figsize=(10, 16), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.0, "wspace": 0.10})
    fig.patch.set_facecolor('xkcd:black')
    new_ax, old_ax = axs.flatten()[0], axs.flatten()[1]
    new_ax.set_title("New EoS")
    old_ax.set_title("Old EoS")
    new_ax.set_ylabel("Temperature (K)")
    for ax in [new_ax, old_ax]:
        ax.set_xlabel(r'Radius $R_{\bigoplus}$')
        ax.grid(alpha=0.4)
    for label in labels.keys():
        new_ax.scatter(
            [p.distance / earth_radius for p in new_particles if p.label == "PLANET" and p.tag == label],
            [p.temperature for p in new_particles if p.label == "PLANET" and p.tag == label],
            label=labels[label],
            s=1
        )
        old_ax.scatter(
            [p.distance / earth_radius for p in old_particles if p.label == "PLANET" and p.tag == label],
            [p.temperature for p in old_particles if p.label == "PLANET" and p.tag == label],
            label=labels[label],
            s=1
        )
    new_mean_curve_x, new_mean_curve_y = __mean_curve(
        x=[p.distance / earth_radius for p in new_particles if p.label == "PLANET"],
        y=[p.temperature for p in new_particles if p.label == "PLANET"],
    )
    old_mean_curve_x, old_mean_curve_y = __mean_curve(
        x=[p.distance / earth_radius for p in old_particles if p.label == "PLANET"],
        y=[p.temperature for p in old_particles if p.label == "PLANET"],
    )
    new_ax.plot(
        new_mean_curve_x, new_mean_curve_y, linewidth=2.0, color='white', label="Mean"
    )
    old_ax.plot(
        old_mean_curve_x, old_mean_curve_y, lioldidth=2.0, color='white', label="Mean"
    )
    legend = new_ax.legend()
    for handle in legend.legendHandles:
        try:
            handle.set_sizes([3.0])
        except:
            pass

    plt.savefig("thermal_profile.png", format='png')


new_particles, new_time = get_particles(path=new_path, number_processes=new_processes, time=time,
                                                solve=False, formatted=True)
old_particles, old_time = get_particles(path=old_path, number_processes=old_processes, time=time,
                                                solve=False, formatted=True)

thermal_profile(new_particles=new_particles, old_particles=old_particles)




