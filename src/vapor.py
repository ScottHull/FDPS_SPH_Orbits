import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

from src.interpolation import NearestNeighbor1D



# def get_particle_vapor_fraction(particle, phase_path):
#     phase_df = pd.read_fwf(phase_path, skiprows=1,
#                            names=["temperature", "density_sol_liq", "density_vap", "pressure",
#                                   "entropy_sol_liq", "entropy_vap"])
#
#     nearest_neighbor = NearestNeighbor1D()
#     nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=particle.temperature,
#                                                                 array=list(phase_df['temperature']))
#     entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
#     entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
#     supercritical = max(phase_df['temperature'])
#     if particle.entropy > supercritical:
#         return 1.0
#     elif particle.entropy < entropy_liq:
#         return 0.0
#     elif entropy_liq <= particle.entropy <= entropy_vap:
#         return (particle.entropy - entropy_liq) / (entropy_vap - entropy_liq)
#     elif particle.entropy > entropy_vap:
#         return 1.0
#
# def get_all_particle_vapor_fractions_from_formatted(df, phase_path, end_label=False):
#     """
#     Returns a list of disk silicate particles with S, T, and their individual VMFs.
#     :param df:
#     :param phase_path:
#     :return:
#     """
#     phase_df = pd.read_fwf(phase_path, skiprows=1,
#                            names=["temperature", "density_sol_liq", "density_vap", "pressure",
#                                   "entropy_sol_liq", "entropy_vap"])
#
#     nearest_neighbor = NearestNeighbor1D()
#     particles = []
#     if end_label:
#         disk_particles = df[df['end_label'] == "DISK"]
#     else:
#         disk_particles = df[df['label'] == "DISK"]
#     sil_disk_particles = disk_particles[disk_particles['tag'] % 2 == 0]
#     s_t = zip(sil_disk_particles['entropy'], sil_disk_particles['temperature'], sil_disk_particles['density'],
#               sil_disk_particles['x'], sil_disk_particles['y'], sil_disk_particles['z'])
#     supercritical = max(phase_df['temperature'])
#     for s, t, rho, x, y, z in s_t:
#         nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
#                                                                     array=list(phase_df['temperature']))
#         entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
#         entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
#         if s > supercritical:
#             return 1.0
#         elif s < entropy_liq:
#             f = 0.0
#         elif entropy_liq <= s <= entropy_vap:
#             f = (s - entropy_liq) / (entropy_vap - entropy_liq)
#         elif s > entropy_vap:
#             f = 1.0
#         particles.append([s, t, f, rho, x, y, z])
#     return particles


def calc_vapor_mass_fraction_from_formatted(df, phase_path, end_label=False):
    phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
    if end_label:
        disk_particles = df[df['end_label'] == "DISK"]
    else:
        disk_particles = df[df['label'] == "DISK"]
    sil_disk_particles = disk_particles[disk_particles['tag'] % 2 == 0]
    s_t = zip(sil_disk_particles['entropy'], sil_disk_particles['temperature'])
    supercritical = max(phase_df['temperature'])
    nearest_neighbor = NearestNeighbor1D()
    num_particles = 0
    vmf = 0
    for s, t in s_t:
        num_particles += 1
        entropy_i = s
        temperature_i = t
        nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                    array=list(phase_df['temperature']))
        entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
        entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
        if t >= supercritical:
            vmf += 1.0
            # particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
        elif s < entropy_liq:
            vmf += 0.0
            # particles['liq'].append((s, t, entropy_liq, entropy_vap))
        elif entropy_liq <= s <= entropy_vap:
            vmf += (s - entropy_liq) / (entropy_vap - entropy_liq)
            # particles['mix'].append((s, t, entropy_liq, entropy_vap))
        elif s > entropy_vap:
            vmf += 1.0
            # particles['vap'].append((s, t, entropy_liq, entropy_vap))
        num_particles += 1

    try:
        vmf = vmf / num_particles
    except Exception as e:  # likely if there are no disk particles
        return 0.0
    return vmf


def calc_vapor_mass_fraction(particles, phase_path, only_disk=True):
    phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    nearest_neighbor = NearestNeighbor1D()
    num_particles = 0
    vmf = 0
    supercritical = max(phase_df['temperature'])
    for p in particles:
        if (p.label == "DISK" or not only_disk) and p.tag % 2 == 0:
            num_particles += 1
            entropy_i = p.entropy
            temperature_i = p.temperature
            nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=temperature_i,
                                                                        array=list(phase_df['temperature']))
            nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=temperature_i,
                                                                        array=list(phase_df['temperature']))
            entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
            entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
            if temperature_i >= supercritical:
                vmf += 1.0
                # particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
            elif entropy_i < entropy_liq:
                vmf += 0.0
                # particles['liq'].append((s, t, entropy_liq, entropy_vap))
            elif entropy_liq <= entropy_i <= entropy_vap:
                vmf += (entropy_i - entropy_liq) / (entropy_vap - entropy_liq)
                # particles['mix'].append((s, t, entropy_liq, entropy_vap))
            elif entropy_i > entropy_vap:
                vmf += 1.0
                # particles['vap'].append((s, t, entropy_liq, entropy_vap))
            num_particles += 1

    try:
        vmf = vmf / num_particles
    except Exception as e:  # likely if there are no disk particles
        return 0.0
    return vmf

def calc_vapor_mass_fraction_with_circularization(particles, phase_path, only_disk=True):
    phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    nearest_neighbor = NearestNeighbor1D()
    num_particles = 0
    vmf = 0
    supercritical = max(phase_df['temperature'])
    for p in particles:
        if (p.label == "DISK" or not only_disk) and p.tag % 2 == 0:
            num_particles += 1
            entropy_i = (p.entropy + p.circularization_entropy_delta)
            temperature_i = p.temperature
            nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=temperature_i,
                                                                        array=list(phase_df['temperature']))
            entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
            entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
            if temperature_i >= supercritical:
                vmf += 1.0
                # particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
            elif entropy_i < entropy_liq:
                vmf += 0.0
                # particles['liq'].append((s, t, entropy_liq, entropy_vap))
            elif entropy_liq <= entropy_i <= entropy_vap:
                vmf += (entropy_i - entropy_liq) / (entropy_vap - entropy_liq)
                # particles['mix'].append((s, t, entropy_liq, entropy_vap))
            elif entropy_i > entropy_vap:
                vmf += 1.0
                # particles['vap'].append((s, t, entropy_liq, entropy_vap))
            num_particles += 1

    try:
        vmf = vmf / num_particles
    except Exception as e:  # likely if there are no disk particles
        return 0.0
    return vmf

def calc_vapor_mass_fraction_with_circularization_from_formatted(df, phase_path, restrict_df=True):
    phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    if restrict_df:
        disk_df = df[df['label'] == "DISK"]
        disk_df = disk_df[disk_df['tag'] % 2 == 0]
        disk_df = disk_df[disk_df['circ_entropy_delta'] < 5000]
        disk_df['entropy_w_circ'] = disk_df['entropy'] + disk_df['circ_entropy_delta']
    else:
        disk_df = df

    nearest_neighbor = NearestNeighbor1D()
    supercritical = max(phase_df['temperature'])
    # get the location of the phase curve that corresponds to the supercritical point
    vmf = 0.0
    num_particles = 0
    # particles = {'liq': [], 'vap': [], 'mix': [], 'supercritical': []}
    for s, t in zip(disk_df['entropy_w_circ'], disk_df['temperature']):
        nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                    array=list(phase_df['temperature']))
        entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
        entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
        if t >= supercritical:
            vmf += 1.0
            # particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
        elif s < entropy_liq:
            vmf += 0.0
            # particles['liq'].append((s, t, entropy_liq, entropy_vap))
        elif entropy_liq <= s <= entropy_vap:
            vmf += (s - entropy_liq) / (entropy_vap - entropy_liq)
            # particles['mix'].append((s, t, entropy_liq, entropy_vap))
        elif s > entropy_vap:
            vmf += 1.0
            # particles['vap'].append((s, t, entropy_liq, entropy_vap))
        num_particles += 1
    try:
        return vmf / num_particles
    except ZeroDivisionError:
        return 0.0

def calc_vapor_mass_fraction_without_circularization_from_formatted(df, phase_path, restrict_df=True):
    phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    if restrict_df:
        disk_df = df[df['label'] == "DISK"]
        disk_df = disk_df[disk_df['tag'] % 2 == 0]
        disk_df = disk_df[disk_df['circ_entropy_delta'] < 5000]
        disk_df['entropy_w_circ'] = disk_df['entropy'] + disk_df['circ_entropy_delta']
    else:
        disk_df = df

    nearest_neighbor = NearestNeighbor1D()
    supercritical = max(phase_df['temperature'])
    vmf = 0.0
    num_particles = 0
    # create a vmf_wo_circ column in the df
    df['vmf_wo_circ'] = 0.0
    # particles = {'liq': [], 'vap': [], 'mix': [], 'supercritical': []}
    for index, s, t in zip(disk_df.index, disk_df['entropy'], disk_df['temperature']):
        nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                    array=list(phase_df['temperature']))
        entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
        entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
        if t >= supercritical:
            vmf += 1.0
            df['vmf_wo_circ'].loc[index] = 1.0
            # particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
        elif s < entropy_liq:
            vmf += 0.0
            df['vmf_wo_circ'].loc[index] = 0.0
            # particles['liq'].append((s, t, entropy_liq, entropy_vap))
        elif entropy_liq <= s <= entropy_vap:
            vmf += (s - entropy_liq) / (entropy_vap - entropy_liq)
            df['vmf_wo_circ'].loc[index] = (s - entropy_liq) / (entropy_vap - entropy_liq)
            # particles['mix'].append((s, t, entropy_liq, entropy_vap))
        elif s > entropy_vap:
            vmf += 1.0
            df['vmf_wo_circ'].loc[index] = 1.0
            # particles['vap'].append((s, t, entropy_liq, entropy_vap))
        num_particles += 1
    try:
        return vmf / num_particles
    except ZeroDivisionError:
        return 0.0


def plot_disk_entropy(particles, phase_path="src/phase_data/duniteS_vapour_curve.txt"):
    phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.distance / 1000 for p in particles if p.label == "DISK" and p.tag % 2 == 0],
        [p.entropy for p in particles if p.label == "DISK" and p.tag % 2 == 0],
        color='black',
        marker="+",
        alpha=0.8
    )
    ax.set_xlabel("Radial Distance (km)")
    ax.set_ylabel("Entropy")
    ax.set_title("Disk Entropy")
    ax.grid()
    plt.savefig("disk_entropy_radial_distance.png", format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.distance / 1000 for p in particles if p.label == "PLANET" and p.tag % 2 == 0],
        [p.entropy for p in particles if p.label == "PLANET" and p.tag % 2 == 0],
        color='black',
        marker="+",
        alpha=0.8
    )
    ax.set_xlabel("Radial Distance (km)")
    ax.set_ylabel("Entropy")
    ax.set_title("Protoplanet Mantle Entropy")
    ax.grid()
    plt.savefig("protoplanet_entropy_radial_distance.png", format='png')
