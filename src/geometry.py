import os
from math import atan, pi, sqrt
from statistics import mean
import matplotlib.pyplot as plt

from src.centering import center_of_mass, center_of_mass_from_formatted


def get_radius(x: list):
    return (max(x) - min(x)) / 2.0


def get_impact_geometry(particles):
    target = [p for p in particles if p.tag <= 1]
    impactor = [p for p in particles if p.tag > 1]

    target_com_x, target_com_y, target_com_z = center_of_mass(
        x_coords=[p.position[0] for p in target],
        y_coords=[p.position[1] for p in target],
        z_coords=[p.position[2] for p in target],
        masses=[p.mass for p in target]
    )
    impactor_com_x, impactor_com_y, impactor_com_z = center_of_mass(
        x_coords=[p.position[0] for p in impactor],
        y_coords=[p.position[1] for p in impactor],
        z_coords=[p.position[2] for p in impactor],
        masses=[p.mass for p in impactor]
    )

    x_offset = impactor_com_x - target_com_x
    y_offset = impactor_com_y - target_com_y
    imp_angle = atan(y_offset / x_offset) * (180.0 / pi)

    min_x_tar = min([p.position[0] for p in target])
    max_x_tar = max([p.position[0] for p in target])

    min_y_tar = min([p.position[1] for p in target])
    max_y_tar = max([p.position[1] for p in target])

    min_x_imp = min([p.position[0] for p in impactor])
    max_x_imp = max([p.position[0] for p in impactor])

    min_y_imp = min([p.position[1] for p in impactor])
    max_y_imp = max([p.position[1] for p in impactor])

    return target, impactor, target_com_x, target_com_y, target_com_z, impactor_com_x, impactor_com_y, impactor_com_z, \
           imp_angle, min_x_tar, max_x_tar, min_y_tar, max_y_tar, min_x_imp, max_x_imp, min_y_imp, max_y_imp


def impact_geometry_image(path, time, iteration, tar_df, imp_df, impactor_com_x, impactor_com_y, target_com_x,
                          target_com_y, imp_angle):
    to_path = path + "_tmp_geometry"
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        tar_df['x'],
        tar_df['y'],
        color='blue',
        s=2,
        label="target"
    )
    ax.scatter(
        imp_df['x'],
        imp_df['y'],
        color='red',
        s=2,
        label="impactor"
    )
    ax.plot(
        [target_com_x, impactor_com_x],
        [target_com_y, target_com_y],
        linewidth=2.0,
        color='white'
    )
    ax.plot(
        [impactor_com_x, impactor_com_x],
        [target_com_y, impactor_com_y],
        linewidth=2.0,
        color='white'
    )
    ax.plot(
        [target_com_x, impactor_com_x],
        [target_com_y, impactor_com_y],
        linewidth=2.0,
        color='white',
        label="IMP ANGLE: {}".format(imp_angle)
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y"),
    ax.set_xlim(-2.1e7, 2.1e7)
    ax.set_ylim(-2.1e7, 2.1e7)
    ax.set_title("Time (hrs): {} // Impact Angle: {}".format(time, imp_angle))
    ax.grid()
    ax.legend(loc='upper left')
    plt.savefig(to_path + "/{}.png".format(iteration), format='png')



def get_velocity_profile(particles, target_radius, impactor_radius):
    G = 6.67 * 10 ** -11
    target = [p for p in particles if p.tag <= 1]
    impactor = [p for p in particles if p.tag > 1]
    target_mass = sum([p.mass for p in target])
    impactor_mass = sum([p.mass for p in impactor])
    total_mass = target_mass + impactor_mass

    v_esc = sqrt((2 * G * total_mass) / (target_radius + impactor_radius))

    target_avg_velocity = [
        mean([p.velocity[0] for p in target]) / v_esc,
        mean([p.velocity[1] for p in target]) / v_esc,
        mean([p.velocity[2] for p in target]) / v_esc,
    ]
    impactor_avg_velocity = [
        mean([p.velocity[0] for p in impactor]) / v_esc,
        mean([p.velocity[1] for p in impactor]) / v_esc,
        mean([p.velocity[2] for p in impactor]) / v_esc,
    ]

    return target_avg_velocity, impactor_avg_velocity, v_esc


def get_impact_geometry_from_formatted(df, name, iteration, time):
    target = df[df['tag'] <= 1]
    impactor = df[df['tag'] > 1]

    target_com_x, target_com_y, target_com_z = center_of_mass_from_formatted(target)
    impactor_com_x, impactor_com_y, impactor_com_z = center_of_mass_from_formatted(impactor)

    x_offset = impactor_com_x - target_com_x
    y_offset = impactor_com_y - target_com_y
    imp_angle = atan(y_offset / x_offset) * (180.0 / pi)

    # min_x_tar = min(target['x'])
    # max_x_tar = max(target['x'])
    #
    # min_y_tar = min(target['y'])
    # max_y_tar = max(target['y'])
    #
    # min_x_imp = min(impactor['x'])
    # max_x_imp = max(impactor['x'])
    #
    # min_y_imp = min(impactor['y'])
    # max_y_imp = max(impactor['y'])

    impact_geometry_image(name, time, iteration, target, impactor, impactor_com_x, impactor_com_y, target_com_x,
                          target_com_y, imp_angle)

    return imp_angle


def get_velocity_profile_from_formatted(df):
    G = 6.67 * 10 ** -11
    target = df[df['tag'] <= 1]
    impactor = df[df['tag'] >= 1]
    target_mass = sum(target['mass'])
    impactor_mass = sum(impactor['mass'])
    total_mass = target_mass + impactor_mass
    target_radius = get_radius(target['x'])
    impactor_radius = get_radius(impactor['x'])

    v_esc = sqrt((2 * G * total_mass) / (target_radius + impactor_radius))

    # target_avg_velocity = [
    #     mean(target['vx']) / v_esc,
    #     mean(target['vy']) / v_esc,
    #     mean(target['vz']) / v_esc,
    # ]
    # impactor_avg_velocity = [
    #     mean(impactor['vx']) / v_esc,
    #     mean(impactor['vy']) / v_esc,
    #     mean(impactor['vz']) / v_esc,
    # ]

    return v_esc
