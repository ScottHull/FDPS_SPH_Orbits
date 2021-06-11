from math import pi, sqrt
import statistics
from copy import copy

from src import elements


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
        statistics.mean([p.relative_velocity[0] for p in particles if p.label == "PLANET" and p.tag == 1]),
        statistics.mean([p.relative_velocity[1] for p in particles if p.label == "PLANET" and p.tag == 1]),
        statistics.mean([p.relative_velocity[1] for p in particles if p.label == "PLANET" and p.tag == 1])
    ]


def collect_particles(output, com, mass_protoearth, relative_velocity=False, find_orbital_elements=True):
    print("Collecting particles...")
    particles = []
    target_velocity = [0, 0, 0]
    if relative_velocity:
        target_velocity = calc_target_velocity(
            vx=output[6],
            vy=output[7],
            vz=output[8],
            tags=output[1]
        )
    for row in output.index:
        position = [float(output[3][row]) - com[0], float(output[4][row]) - com[1],
                    float(output[5][row]) - com[2]]
        velocity = [float(output[6][row]), float(output[7][row]), float(output[8][row])]
        relative_velocity = [
            velocity[0] - target_velocity[0],
            velocity[1] - target_velocity[1],
            velocity[2] - target_velocity[2]
        ]
        p = elements.Particle(
            particle_id=int(output[0][row]),
            tag=int(output[1][row]),
            mass=float(output[2][row]),
            position=position,
            velocity=velocity,
            relative_velocity=relative_velocity,
            density=float(output[9][row]),
            internal_energy=float(output[10][row]),
            pressure=float(output[11][row]),
            potential_energy=float(output[12][row]),
            entropy=float(output[13][row]),
            temperature=float(output[14][row]),
            mass_grav_body=float(mass_protoearth),
            calculate_elements=find_orbital_elements
        )
        particles.append(p)
    print("Collected particles!")
    return particles


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


def will_be_planet_circular_orbit(p, a):
    if abs(p.radius_circular_orbit) < a:
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
    if p.eccentricity <= 1.0 and abs(p.periapsis) > a:
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
