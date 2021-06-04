from math import pi, sqrt
import statistics
from copy import copy


def calc_oblateness(a, b):
    return (a - b) / b


def refined_oblateness(T_star, T_protoearth, K=0.335):
    numerator = (5.0 / 2.0) * ((T_star / T_protoearth) ** 2)
    denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
    return numerator / denominator


def calc_mass_protoearth(a, b, rho=5500):
    return ((4.0 / 3.0) * pi * (a ** 2) * b) * rho


def calc_target_velocity(vx, vy, vz, tags):
    return [
        statistics.mean([i for index, i in enumerate(vx) if tags[index] == 1]),
        statistics.mean([i for index, i in enumerate(vy) if tags[index] == 1]),
        statistics.mean([i for index, i in enumerate(vz) if tags[index] == 1])
    ]


def refine_target_velocity(particles):
    # returns relative velocity of target iron
    return [
        statistics.mean([p.relative_velocity[0] for p in particles if particles.label == "PLANET" and p.tag == 1]),
        statistics.mean([p.relative_velocity[1] for p in particles if particles.label == "PLANET" and p.tag == 1]),
        statistics.mean([p.relative_velocity[1] for p in particles if particles.label == "PLANET" and p.tag == 1])
    ]


def is_planet(p, a):
    """
    If the particle falls within the radius of the planet.
    :param p:
    :param a:
    :return:
    """
    if p.distance < a:
        p.label = "PLANET"
        return True
    return False


def will_be_planet(p, a):
    """
    If the particle's periapsis is within the radius of the planet.
    :param p:
    :param a:
    :return:
    """
    if p.eccentricity < 1.0 and abs(p.periapsis) <= a:
        p.label = "PLANET"
        return True
    return False


def is_disk(p, a):
    if p.eccentricity < 1.0 and abs(p.periapsis) > a:
        p.label = "DISK"
        return True
    return False


def is_escape(p, a):
    if p.eccentricity > 1.0:
        p.label = "ESCAPE"
        return True
    return False


def log(iteration, error, a,
        NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
        NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING, NEW_MASS_PROTOPLANET, NEW_MASS_DISK, NEW_MASS_ESCAPED,
        total_angular_momentum):
    EARTH_MASS = 5.972 * 10 ** 24
    LUNAR_MASS = 7.34767309 * 10 ** 22
    L_EM = 3.5 * 10 ** 34
    print(
        "ITERATION: {}\n"
        "ERROR: {}\n"
        "NEW A: {}\n"
        "NUM_PARTICLES_WITHIN_RADIAL_DISTANCE: {}\n"
        "NUM_PARTICLES_IN_DISK: {}\n"
        "NUM_PARTICLES_ESCAPING: {}".format(iteration, error, a,
                                            NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
                                            NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING)
    )
    print(
        "PROTOPLANET MASS: {} M_E ({} KG)\n"
        "DISK MASS: {} M_L ({} KG)\n"
        "ESCAPING MASS: {} M_L ({} KG)\n".format(
            NEW_MASS_PROTOPLANET / EARTH_MASS, NEW_MASS_PROTOPLANET,
            NEW_MASS_DISK / LUNAR_MASS, NEW_MASS_DISK,
            NEW_MASS_ESCAPED / LUNAR_MASS, NEW_MASS_ESCAPED
        )
    )
    print(
        "TOTAL ANGULAR MOMENTUM: {} L_EM ({})\n\n".format(total_angular_momentum / L_EM, total_angular_momentum)
    )
