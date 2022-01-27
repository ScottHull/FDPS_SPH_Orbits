import os
import csv
import shutil
from math import pi, asin, isnan
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile
from src.theia import LunaToTheia


def get_time(f, local=True):
    formatted_time = None
    if local:  # if not reading from remote server
        with open(f, 'r') as infile:
            reader = csv.reader(infile, delimiter="\t")
            formatted_time = float(next(reader)[0])
        infile.close()
    else:
        formatted_time = float(next(f)[0])
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours


def get_setup_file_data(path):
    d = {}
    try:
        with open(path, 'r') as infile:
            found_initial_setup = False
            for line in infile:
                if not found_initial_setup:
                    if "INITIAL SETUP" in line:
                        found_initial_setup = True
                else:
                    if len(line) == 1:
                        break
                    else:
                        header, data = line.split(":")
                        d.update({header: data.replace("\n", "")})
    except:
        pass
    return d


def __get_vmf_timeplot_data(path, phase_path, start_iteration, end_iteration, increment, L_EM=3.5 * 10 ** 34):
    max_time = get_time(path + "/{}.csv".format(end_iteration))
    times = []
    vmfs = []
    disk_particle_count = []
    avg_disk_entropy = []
    spec_disk_ams = []
    for time in np.arange(start_iteration, end_iteration + increment, increment):
        f = path + "/{}.csv".format(time)
        try:
            times.append(get_time(f))
            df = pd.read_csv(f, skiprows=2)
            vmf = calc_vapor_mass_fraction_from_formatted(df=df, phase_path=phase_path) * 100.0
            vmfs.append(vmf)

            disk_particles = df[df['label'] == "DISK"]
            positions = list(zip(disk_particles['x'], disk_particles['y'], disk_particles['z']))
            velocities = list(zip(disk_particles['vx'], disk_particles['vy'], disk_particles['vz']))
            masses = list(disk_particles['mass'])
            spec_disk_ams.append(sum([masses[index] * np.linalg.norm(np.cross(p, velocities[index])) for index, p in
                                      enumerate(positions)]) / L_EM)  # specific angular momentum of the disk
            try:
                avg_disk_entropy_at_time = mean(disk_particles['entropy'])
            except:
                avg_disk_entropy_at_time = 0
            avg_disk_entropy.append(avg_disk_entropy_at_time)
            num_disk_particles = len(disk_particles['entropy'])
            disk_particle_count.append(num_disk_particles)
        except:
            times.append(np.nan)
            vmfs.append(np.nan)
            disk_particle_count.append(np.nan)
            avg_disk_entropy.append(np.nan)
    return times, vmfs, disk_particle_count, avg_disk_entropy, max_time, spec_disk_ams


def build_vmf_timeplots(name, meta, start_iteration, end_iteration, increment, label_header='label',
                        output_fname="vmf_profile_output.txt", L_EM=3.5 * 10 ** 34):
    """
    Builds VMF timeseries plots from formatted outputs.
    :param names:
    :param paths:
    :return:
    """
    plt.style.use("dark_background")
    fig, axs = plt.subplots(4, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    axs = axs.flatten()
    axs[0].set_title("New EoS")
    axs[1].set_title("Old EoS")
    for ax in axs:
        ax.grid(alpha=0.4)
    for i in meta.keys():
        try:
            print("working on {}".format(i))
            n = meta[i]['name']
            p = meta[i]['path']
            if "new" in p:
                phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
            else:
                phase_path = "src/phase_data/duniteN__vapour_curve.txt"
            times, vmfs, disk_particle_count, avg_disk_entropy, max_time, spec_ang_mom = __get_vmf_timeplot_data(
                p, phase_path, start_iteration, end_iteration, increment, L_EM=L_EM)
            if "new" in i:
                axs[0].plot(times, vmfs, linewidth=2.0, label=n)
                axs[0].set_ylabel("VMF (%)")
                axs[2].plot(times, avg_disk_entropy, linewidth=2.0, label=n)
                axs[2].set_ylabel("Avg. Disk Entropy")
                axs[4].plot(times, disk_particle_count, linewidth=2.0, label=n)
                axs[4].set_ylabel("# Disk Particles")
                axs[6].plot(times, spec_ang_mom, linewidth=2.0, label=n)
                axs[6].set_ylabel("Disk Angular Momentum (L_EM)")
            else:
                axs[1].plot(times, vmfs, linewidth=2.0, label=n)
                axs[3].plot(times, avg_disk_entropy, linewidth=2.0, label=n)
                axs[5].plot(times, disk_particle_count, linewidth=2.0, label=n)
                axs[7].plot(times, spec_ang_mom, linewidth=2.0, label=n)
        except FileNotFoundError:
            print(i)
    for ax in axs:
        ax.legend(loc='lower right')
    axs[0].set_ylim(0, 45)
    axs[1].set_ylim(0, 45)
    axs[2].set_ylim(0, 7000)
    axs[3].set_ylim(0, 7000)
    axs[4].set_ylim(0, 35000)
    axs[5].set_ylim(0, 35000)
    axs[6].set_ylim(0, 0.5)
    axs[7].set_ylim(0, 0.5)
    axs[-2].set_xlabel("Time (hrs)")
    axs[-1].set_xlabel("Time (hrs)")
    plt.savefig("vmf_timeseries_{}.png".format(name), format='png', dpi=200)


def build_impact_angle_geometries(name, meta, start_iteration, end_iteration, specified_imp_angle, increment=1):
    plt.style.use("dark_background")
    """
    :param names:
    :param paths:
    :param start_iteration:
    :param end_iteration:
    :param increment:
    :return:
    """

    d = {}

    for i in meta.keys():
        try:
            print("working on {}".format(i))
            n = meta[i]['name']
            p = meta[i]['path']
            times, imp_angles = [], []
            for time in np.arange(start_iteration, end_iteration + increment, increment):
                t = get_time(p + "/{}.csv".format(time))
                times.append(t)
                df = pd.read_csv(p + "/{}.csv".format(time), skiprows=2)
                imp_angles.append(get_impact_geometry_from_formatted(df=df, name=i, time=t, iteration=time))
            d.update({n: {"imp_angles": imp_angles, "times": times}})
            animate(
                start_time=start_iteration,
                end_time=end_iteration,
                interval=increment,
                path=i + "_tmp_geometry",
                filename="{}_imp_geometry.mp4".format(i)
            )
        except FileNotFoundError:
            print(i)

    imp_ang_fig, imp_ang_axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    imp_ang_axs = imp_ang_axs.flatten()
    imp_ang_axs[0].set_title("Impact Parameter (New EoS)"), imp_ang_axs[1].set_title("Impact Parameter (Old EoS)")
    imp_ang_axs[0].set_ylabel("Impact Angle (deg)")

    for ax in imp_ang_axs:
        ax.axhline(asin(specified_imp_angle) * (180 / pi), linewidth=2.0, linestyle="--", c='red',
                   label="{} deg.".format(round(asin(specified_imp_angle) * (180 / pi), 2)))
        ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)

    for n in d.keys():
        if "n" in n.lower():
            imp_ang_axs[0].plot(
                d[n]["times"], d[n]["imp_angles"], linewidth=2.0, label=n
            )
        else:
            imp_ang_axs[1].plot(
                d[n]["times"], d[n]["imp_angles"], linewidth=2.0, label=n
            )
    for ax in imp_ang_axs:
        ax.legend(loc='upper left')
    plt.savefig("impact_angle_profile_{}.png".format(name), format='png', dpi=200)


def build_impact_velocity_charts(name, meta, start_iteration, end_iteration, increment=1):
    plt.style.use("dark_background")
    imp_vel_fig, imp_vel_axs = plt.subplots(len(meta.keys()), 1, figsize=(16, 32), sharex='all', sharey='all',
                                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    imp_vel_axs = imp_vel_axs.flatten()

    fig_index = 0

    for ax in imp_vel_axs:
        ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)

    for i in meta.keys():
        try:
            print("working on {}".format(i))
            n = meta[i]['name']
            p = meta[i]['path']
            specified_imp_vel = float(meta[i]['setup']["ESCAPE VELOCITY"])
            times, imp_vels = [], []
            for time in np.arange(start_iteration, end_iteration + increment, increment):
                times.append(get_time(p + "/{}.csv".format(time)))
                df = pd.read_csv(p + "/{}.csv".format(time), skiprows=2)
                imp_vels.append(get_velocity_profile_from_formatted(df) / 1000)
            imp_vel_axs[fig_index].plot(
                times, imp_vels, linewidth=2.0,
                label="{} (Max: {})".format(n, round(max([x for x in imp_vels if isnan(x) == False]), 3))
            )
            imp_vel_axs[fig_index].axhline(specified_imp_vel / 1000, color='red', linewidth=2.0, linestyle="--",
                                           label="Specified ({} km/s)".format(round(specified_imp_vel / 1000, 3)))
            imp_vel_axs[fig_index].set_ylabel("Impact Velocity (km/s)")
            fig_index += 1
        except FileNotFoundError:
            print(i)
    for ax in imp_vel_axs:
        ax.legend(loc='upper left')
    plt.savefig("impact_velocity_profile_{}.png".format(name), format='png', dpi=200)


def map_disk_to_phase_profile(name, meta, end_iteration):
    plt.style.use("dark_background")
    new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
    old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
    new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    fig, axs = plt.subplots(len(meta.keys()), 1, figsize=(16, 32), sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    axs = axs.flatten()
    fig_index = 0
    for ax in axs:
        ax.set_xlim(0, 15000)
        # ax.set_ylim(0, 15000)
        ax.grid(alpha=0.4)

    for i in meta.keys():
        try:
            n = meta[i]['name']
            p = meta[i]['path']
            df = pd.read_csv(p + "/{}.csv".format(end_iteration), skiprows=2)
            disk = df[df['tag'] % 2 == 0]
            disk = disk[disk['label'] == "DISK"]
            if "new" in i.lower():
                axs[fig_index].plot(
                    new_phase_df['entropy_vap'],
                    new_phase_df['temperature'],
                    linewidth=2.0,
                    label="Vapor"
                )
                axs[fig_index].plot(
                    new_phase_df['entropy_sol_liq'],
                    new_phase_df['temperature'],
                    linewidth=2.0,
                    label="Liquid",
                )
            else:
                axs[fig_index].plot(
                    old_phase_df['entropy_vap'],
                    old_phase_df['temperature'],
                    linewidth=2.0,
                    label="Vapor"
                )
                axs[fig_index].plot(
                    old_phase_df['entropy_sol_liq'],
                    old_phase_df['temperature'],
                    linewidth=2.0,
                    label="Liquid",
                )
            axs[fig_index].scatter(
                disk['entropy'],
                disk['temperature'],
                s=2,
                c='#fa8174',
                label="{} disk particles".format(n)
            )
            axs[fig_index].set_ylabel("Temperature")
            axs[fig_index].legend(loc='upper left')
            fig_index += 1
        except Exception as e:
            print("Problem!", e)
    axs[-1].set_xlabel("Entropy")
    axs[0].set_title("Disk Particles on Phase Curve")
    plt.savefig("disk_on_phase_curve_{}.png".format(name), format='png', dpi=200)


def map_disk_to_phase_profile_eos_charts(name, meta, end_iteration):
    plt.style.use("dark_background")
    new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
    old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
    new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    axs = axs.flatten()
    axs[0].set_title("New EoS")
    axs[1].set_title("Old EoS")
    axs[0].plot(
        new_phase_df['entropy_vap'],
        new_phase_df['temperature'],
        linewidth=2.0,
        label="Vapor"
    )
    axs[0].plot(
        new_phase_df['entropy_sol_liq'],
        new_phase_df['temperature'],
        linewidth=2.0,
        label="Liquid",
    )
    axs[1].plot(
        old_phase_df['entropy_vap'],
        old_phase_df['temperature'],
        linewidth=2.0,
        label="Vapor"
    )
    axs[1].plot(
        old_phase_df['entropy_sol_liq'],
        old_phase_df['temperature'],
        linewidth=2.0,
        label="Liquid",
    )
    axs[0].set_ylabel("Temperature")
    for ax in axs:
        ax.set_xlim(0, 15000)
        # ax.set_ylim(0, 15000)
        ax.grid(alpha=0.4)
        ax.set_xlabel("Entropy")

    for i in meta.keys():
        try:
            fig_index = None
            n = meta[i]['name']
            p = meta[i]['path']
            df = pd.read_csv(p + "/{}.csv".format(end_iteration), skiprows=2)
            disk = df[df['tag'] % 2 == 0]
            disk = disk[disk['label'] == "DISK"]
            if "new" in i.lower():
                fig_index = 0
            else:
                fig_index = 1
            axs[fig_index].scatter(
                disk['entropy'],
                disk['temperature'],
                s=2,
                label="{} disk particles".format(n)
            )
        except Exception as e:
            print("Problem!", e)
    for ax in axs:
        ax.legend(loc='upper left')
    plt.savefig("disk_on_phase_curve_same_eos_plots_{}.png".format(name), format='png', dpi=200)


def __profile_time(a):
    meta, i, end_iteration, number_processes = a
    new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
    old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
    try:
        if "new" in i:
            phase_path = new_phase_path
        else:
            phase_path = old_phase_path
        fig_index = None
        n = meta[i]['name']
        p = meta[i]['path'].replace("formatted_", "")
        to_fname = "merged_{}_{}.dat".format(end_iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=end_iteration, output_path=p, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = cf.sim_time
        f = os.getcwd() + "/{}".format(to_fname)
        pm = ParticleMap(path=f, center=True, relative_velocity=False)
        particles = pm.collect_particles(find_orbital_elements=True)
        pm.solve(particles=particles, phase_path=phase_path, report_name="{}.txt".format(i),
                 iteration=end_iteration, simulation_time=formatted_time)
        os.remove(to_fname)
    except Exception as e:
        print("problem!", e)
    return None


def get_end_profile_reports(meta, end_iteration, number_processes=200):
    pool = mp.Pool(5)
    pool.map(__profile_time, [[meta, i, end_iteration, number_processes] for i in meta.keys()])
    pool.close()
    pool.join()


def __build_scene(d):
    meta, iteration, to_path, min_normalize_parameter, max_normalize_parameter, square_scale, s, u, p, to_client_path = d
    client = LunaToTheia(s, u, p)
    normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
    cmap = cm.get_cmap('jet')

    num_new = len([i for i in meta.keys() if "new" in i])
    num_old = len([i for i in meta.keys() if "old" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-square_scale, square_scale)
        ax.set_ylim(-square_scale, square_scale)
        ax.set_xticks([], minor=False)
        ax.set_yticks([], minor=False)
    index_new, index_old = 0, 1
    for i in meta.keys():
        n = meta[i]['name']
        p = meta[i]['path']
        if s is not None:
            f = client.get_file(client.theia_client, p, "{}.csv".format(iteration))
            formatted_time = get_time(f, local=False)
            df = client.get_df_from_theia(p, "{}.csv".format(iteration), skiprows=2)
        else:
            formatted_time = get_time(p + "/{}.csv".format(iteration))
            df = pd.read_csv(p + "/{}.csv".format(iteration), skiprows=2)
        df = df[df['z'] < 0]
        if "new" in i:
            axs[index_new].scatter(
                df['x'], df['y'], s=1, color=[cmap(normalizer(i)) for i in df['entropy']]
            )
            axs[index_new].set_title(n + " {} hrs".format(formatted_time))
            index_new += 2
        else:
            axs[index_old].scatter(
                df['x'], df['y'], s=1, color=[cmap(normalizer(i)) for i in df['entropy']]
            )
            axs[index_old].set_title(n + " ({} hrs)".format(formatted_time))
            index_old += 2
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(axs[0], width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_title("Entropy", fontsize=6)
    plt.savefig(to_path + "/{}.png".format(iteration), format='png')
    if s is not None:
        client.send_file_to_theia(to_path, to_client_path, "/{}.png".format(iteration))
        os.remove(to_path + "/{}.png".format(iteration))
    client.theia_client.close()


def build_scenes(name, meta, to_path, start_iteration, end_iteration, increment, s=None, u=None, p=None, to_client_path="", fill=False, proc=10):
    # if os.path.exists(to_path):
    #     shutil.rmtree(to_path)
    if not fill:
        try:
            os.mkdir(to_path)
        except:
            pass
    pool = mp.Pool(proc)
    to_make = np.arange(start_iteration, end_iteration + increment, increment)
    if fill:
        if s is not None:
            client = LunaToTheia(s, u, p)
            ldir = client.listdir(client.theia_client, to_client_path)
            to_make = [i for i in to_make if "{}.png".format(i) not in ldir]
        else:
            to_make = [i for i in to_make if "{}.png".format(i) not in os.listdir(to_path)]
        to_make = [i for i in to_make if "{}.png".format(i) not in os.listdir(to_path)]
    pool.map(__build_scene, [[meta, iteration, to_path, 2000, 8000, 4e7, s, u, p, to_client_path] for iteration in
                             to_make])
    pool.close()
    pool.join()
    animate(
        start_time=start_iteration,
        end_time=end_iteration,
        interval=increment,
        path=to_path,
        fps=30,
        filename="{}_animated_from_formatted.mp4".format(name),
    )


def disk_temperature_vs_radius(name, meta, iteration):
    plt.style.use("dark_background")
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    axs = axs.flatten()
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.set_xlabel(r"Distance from Earth Center (Radius $R_{\bigoplus}$)")
    axs[0].set_ylabel("Temperature (K)")
    fig.patch.set_facecolor('xkcd:black')
    for i in meta.keys():
        try:
            n = meta[i]['name']
            p = meta[i]['path']
            axs[0].set_title("New EoS ({} hrs)".format(get_time(p + "/{}.csv".format(iteration))))
            axs[1].set_title("Old EoS ({} hrs)".format(get_time(p + "/{}.csv".format(iteration))))
            df = pd.read_csv(p + "/{}.csv".format(iteration), skiprows=2)
            disk = df[df['tag'] % 2 == 0]
            disk = disk[disk['label'] == "DISK"]
            if "new" in i:
                axs[0].scatter(
                    disk['radius'] / (6371 * 1000),
                    disk['temperature'],
                    s=2,
                    label=n
                )
            else:
                axs[1].scatter(
                    disk['radius'] / (6371 * 1000),
                    disk['temperature'],
                    s=2,
                    label=n
                )
        except Exception as e:
            print(e)
            pass
    for ax in axs:
        ax.legend(loc='upper right')
    fig.suptitle(name)

    plt.savefig("{}_disk_temperatures.png".format(name), format='png', dpi=200)
