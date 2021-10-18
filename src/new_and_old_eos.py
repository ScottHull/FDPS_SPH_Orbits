#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times


def get_particles(path, number_processes, time, solve=False):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    # f = os.getcwd() + "/merged_{}.dat".format(closest_iteration_to_time)
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    if solve:
        pm.solve(particles=particles)
    os.remove(f)
    return particles, formatted_time


def get_parameter_from_particles(particle, parameter):
    d = {
        "entropy": particle.entropy,
        "internal_energy": particle.internal_energy,
        "temperature": particle.temperature,
    }
    return d[parameter]


def plot(fig, axs, particles, index, time, cmap, normalizer, square_scale, parameter, flatten=True):
    if flatten:
        ax = axs.flatten()[index]
    else:
        ax = axs[index]
    # ax.set_facecolor('xkcd:black')
    # ax.spines['left'].set_color('white')
    # ax.spines['right'].set_color('white')
    # ax.spines['bottom'].set_color('white')
    # ax.spines['top'].set_color('white')
    ax.scatter(
        [p.position[0] for p in particles if p.position[2] < 0],
        [p.position[1] for p in particles if p.position[2] < 0],
        s=0.02,
        marker="o",
        c=[cmap(normalizer(get_parameter_from_particles(particle=p, parameter=parameter))) for p in particles if
           p.position[2] < 0],
    )
    ax.text(
        square_scale - (square_scale / 1.2),
        square_scale - (square_scale / 3),
        str(round(seconds_to_hours(time), 2)) + " hrs",
        c="white",
        fontsize=10
    )
    # ax.set_title(seconds_to_hours(time))
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)
    ax.set_yticks([])
    # for minor ticks
    ax.set_yticks([], minor=True)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_box_aspect(1)

    scalebar = AnchoredSizeBar(ax.transData,
                               6378.1 * 1000,
                               r'1 $R_{\bigoplus}$',
                               loc=8,
                               pad=0.3,
                               color='white',
                               frameon=False,
                               size_vertical=1,
                               fontproperties=fm.FontProperties(size=6))

    ax.add_artist(scalebar)
    return ax


def scatter(fig, axs, particles, index, square_scale, parameter):
    ax = axs.flatten()[index]
    ax.scatter(
        [p.distance for p in particles if p.position[2] < 0 and p.tag == 0],
        [get_parameter_from_particles(particle=p, parameter=parameter) for p in particles if
         p.position[2] < 0 and p.tag == 0],
        s=0.02,
        marker="o",
        label="Target Silicate"
    )
    ax.scatter(
        [p.distance for p in particles if p.position[2] < 0 and p.tag == 1],
        [get_parameter_from_particles(particle=p, parameter=parameter) for p in particles if
         p.position[2] < 0 and p.tag == 1],
        s=0.02,
        marker="o",
        label="Target Iron"
    )
    ax.scatter(
        [p.distance for p in particles if p.position[2] < 0 and p.tag == 2],
        [get_parameter_from_particles(particle=p, parameter=parameter) for p in particles if
         p.position[2] < 0 and p.tag == 2],
        s=0.02,
        marker="o",
        label="Impactor Silicate"
    )
    ax.scatter(
        [p.distance for p in particles if p.position[2] < 0 and p.tag == 3],
        [get_parameter_from_particles(particle=p, parameter=parameter) for p in particles if
         p.position[2] < 0 and p.tag == 3],
        s=0.02,
        marker="o",
        label="Impactor Iron"
    )
    ax.set_xlim(0, square_scale)
    ax.set_box_aspect(1)
    return ax


def main_plotting_loop(min_iteration, max_iteration, number_processes, time, new_path, old_path, normalizer,
                       square_scale, cmap, to_path, parameter):
    new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=time)
    old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=time)
    fig, axs = plt.subplots(2, 2, figsize=(10, 10),
                            gridspec_kw={"hspace": 0.0, "wspace": 0.08})
    fig.patch.set_facecolor('xkcd:black')
    ax1 = plot(fig=fig, axs=axs, index=0, time=new_time, particles=new_particles, cmap=cmap,
               normalizer=normalizer,
               parameter=parameter, square_scale=square_scale)
    ax2 = plot(fig=fig, axs=axs, index=1, time=new_time, particles=old_particles, cmap=cmap,
               normalizer=normalizer,
               parameter=parameter, square_scale=square_scale)
    ax3 = scatter(fig=fig, axs=axs, index=2, particles=new_particles,
                  parameter=parameter, square_scale=square_scale)
    ax4 = scatter(fig=fig, axs=axs, index=3, particles=old_particles,
                  parameter=parameter, square_scale=square_scale)
    axs.flatten()[2].set_ylabel(parameter.replace("_", " ").title())

    axs.flatten()[0].set_title("New EoS")
    axs.flatten()[1].set_title("Old EoS")
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    # cbar = fig.colorbar(sm, ax=axs.flatten()[1])
    cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    # cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.set_title(parameter.replace("_", " ").title(), fontsize=6)
    legend = ax3.legend(fontsize=6, loc='upper left')
    for handle in legend.legendHandles:
        handle.set_sizes([3.0])
    ax4.set_yticks([])
    ax4.set_yticks([], minor=True)
    ax3_ymin, ax3_ymax = ax3.get_ylim()
    ax4_ymin, ax4_ymax = ax4.get_ylim()
    lims = [ax3_ymin, ax3_ymax, ax4_ymin, ax4_ymax]
    scatter_range = [min(lims), max(lims)]
    inc = (scatter_range[1] - scatter_range[0]) * 0.1
    ax3.set_ylim(0, scatter_range[1] + inc)
    ax4.set_ylim(0, scatter_range[1] + inc)
    ax3.set_xlabel("Radius (m)")
    ax4.set_xlabel("Radius (m)")
    plt.savefig(to_path + "/{}.png".format(time), format='png', dpi=200)


def get_parameter(particle, time, iteration):
    return {
        "pressure": particle.pressure,
        "internal_energy": particle.internal_energy,
        "entropy": particle.entropy,
        "temperature": particle.temperature,
        "density": particle.density,
        "tag": particle.tag,
        "time": time,
        "iteration": iteration
    }

def __peaks_df(peaks, name, save=False):
    particle_ids = list(peaks.keys())
    df = pd.DataFrame({
        "particle_id": particle_ids,
        "pressure": [peaks[p_id]["pressure"] for p_id in particle_ids],
        "internal_energy": [peaks[p_id]["internal_energy"] for p_id in particle_ids],
        "entropy": [peaks[p_id]["entropy"] for p_id in particle_ids],
        "temperature": [peaks[p_id]["temperature"] for p_id in particle_ids],
        "density": [peaks[p_id]["density"] for p_id in particle_ids],
        "tag": [peaks[p_id]["tag"] for p_id in particle_ids],
        "time": [peaks[p_id]["time"] for p_id in particle_ids],
        "iteration": [peaks[p_id]["iteration"] for p_id in particle_ids],
    })
    if save:
        df.to_csv("peaks_{}.csv".format(name))
    return df

def get_peak(save, parameter, min_iteration, max_iteration, interval, new_path, old_path, number_processes, solve=False):
    d_new = {}
    d_old = {}
    __iter = 0
    for time in np.arange(min_iteration, max_iteration + interval, interval):
        new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=time,
                                                solve=solve)
        old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=time,
                                                solve=solve)
        if __iter == 0:
            for p in new_particles:
                d_new.update({p.particle_id: get_parameter(particle=p, time=new_time, iteration=time)})
            for p in old_particles:
                d_old.update({p.particle_id: get_parameter(particle=p, time=old_time, iteration=time)})
        for p in new_particles:
            params = get_parameter(particle=p, time=new_time, iteration=time)
            current = d_new[p.particle_id]
            if params[parameter] > current[parameter]:
                d_new[p.particle_id] = params
        for p in old_particles:
            params = get_parameter(particle=p, time=old_time, iteration=time)
            current = d_old[p.particle_id]
            if params[parameter] > current[parameter]:
                d_old[p.particle_id] = params
        __iter += 1
    peaks_new_df = __peaks_df(peaks=d_new, name="new", save=save)
    peaks_old_df = __peaks_df(peaks=d_old, name="old", save=save)
    return d_new, d_old
