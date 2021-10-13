#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times


def get_particles(path, number_processes, time):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    # f = os.getcwd() + "/merged_{}.dat".format(closest_iteration_to_time)
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    # pm.solve(particles=particles)
    os.remove(f)
    return particles, formatted_time


def get_parameter_from_particles(particle, parameter):
    d = {
        "entropy": particle.entropy,
        "internal_energy": particle.internal_energy,
        "temperature": particle.temperature,
    }
    return d[parameter]


def plot(fig, axs, particles, index, time, cmap, normalizer, square_scale, parameter):
    ax = axs.flatten()[index]
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
                               square_scale / 5,
                               '{:.2e} km'.format(square_scale / 5),
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
