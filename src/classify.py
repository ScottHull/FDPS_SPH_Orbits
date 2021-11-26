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


def get_iron_fraction(particles):
    try:
        roche_limit = get_roche_radius()
        disk_particles = [p for p in particles if p.label == "DISK"]
        total_disk_mass = sum([p.mass for p in disk_particles])
        total_iron_disk_mass = sum([p.mass for p in disk_particles if p.tag % 2 != 0])
        total_iron_disk_mass_beyond_roche = sum(
            [p.mass for p in disk_particles if p.tag % 2 != 0 and p.distance > roche_limit])
        return total_iron_disk_mass / total_disk_mass, total_iron_disk_mass_beyond_roche / total_disk_mass
    except Exception as e:
        return 0.0, 0.0


def refine_target_velocity(particles):
    # returns relative velocity of target iron
    return [
        statistics.mean([p.relative_velocity[0] for p in particles if p.label == "PLANET" and p.tag == 1]),
        statistics.mean([p.relative_velocity[1] for p in particles if p.label == "PLANET" and p.tag == 1]),
        statistics.mean([p.relative_velocity[1] for p in particles if p.label == "PLANET" and p.tag == 1])
    ]


def collect_particles(output, com, mass_protoearth, relative_velocity=False, find_orbital_elements=True, formatted=False):
    """
    Columns for formatted files:
    id,x,y,z,vx,vy,vz,mass,radius,tag,label,entropy,temperature,density,pressure,internal_energy
    :param output:
    :param com:
    :param mass_protoearth:
    :param relative_velocity:
    :param find_orbital_elements:
    :param formatted:
    :return:
    """
    output = output.to_dict()
    print("Collecting particles...")
    particles = []
    target_velocity = [0, 0, 0]
    if relative_velocity:
        if not formatted:
            target_velocity = calc_target_velocity(
                vx=output[6],
                vy=output[7],
                vz=output[8],
                tags=output[1]
            )
        else:
            target_velocity = calc_target_velocity(
                vx=output["vx"],
                vy=output["vy"],
                vz=output["vz"],
                tags=output["tag"]
            )
    if not formatted:
        for row in range(0, len(output[list(output.keys())[0]])):
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
    else:
        for row in range(0, len(output[list(output.keys())[0]])):
            position = [float(output["x"][row]) - com[0], float(output["y"][row]) - com[1],
                        float(output["z"][row]) - com[2]]
            velocity = [float(output["vx"][row]), float(output["vy"][row]), float(output["vz"][row])]
            relative_velocity = [
                velocity[0] - target_velocity[0],
                velocity[1] - target_velocity[1],
                velocity[2] - target_velocity[2]
            ]
            p = elements.Particle(
                particle_id=int(output["id"][row]),
                tag=int(output["tag"][row]),
                mass=float(output["mass"][row]),
                position=position,
                velocity=velocity,
                relative_velocity=relative_velocity,
                density=float(output["density"][row]),
                internal_energy=float(output["internal_energy"][row]),
                pressure=float(output["pressure"][row]),
                potential_energy=float(output["potential_energy"][row]),
                entropy=float(output["entropy"][row]),
                temperature=float(output["temperature"][row]),
                mass_grav_body=float(mass_protoearth),
                calculate_elements=find_orbital_elements
            )
            p.label = output['label'][row]
            particles.append(p)
    print("Collected particles!")
    return particles


def average_density(planet_mass, a):
    vol_sphere = (4 / 3) * pi * (a ** 3)
    return planet_mass / vol_sphere


def get_roche_radius():
    radius_earth = 6371 * 1000
    return 2.9 * radius_earth


def predicted_satellite_mass(disk_angular_momentum, mass_target, mass_disk, mass_escape):
    # Canup 2004 equation 1
    G = 6.67 * 10 ** -11
    radius_earth = 6371 * 1000
    mass_earth = 5.972e24
    roche_radius = get_roche_radius()
    mass_escape = 0.05 * mass_disk  # assumption based on Canup 2004 for more centrally condensed disks than Ida 1997
    a1 = 1.9 * disk_angular_momentum / sqrt(G * mass_earth * roche_radius)
    a2 = 1.1 * mass_disk
    a3 = 1.9 * mass_escape
    return a1 - a2 - a3


def is_beyond_roche_radius(p):
    roche_radius = get_roche_radius()
    if p.distance > roche_radius:
        return True
    return False


def is_planet(p, a):
    """
    If the particle falls within the radius of the planet.
    :param p:
    :param a:
    :return:
    """
    try:
        if abs(p.distance) < a:
            p.label = "PLANET"
            return True
        return False
    except Exception as e:
        p.label = "ERROR"
        return False


def circular_orbit_beyond_roche(p):
    roche = get_roche_radius()
    if abs(p.radius_circular_orbit) > roche:
        return True
    return False


def will_be_planet_circular_orbit(p, a):
    try:
        if abs(p.radius_circular_orbit) < a:
            p.label = "PLANET"
            return True
        return False
    except Exception as e:
        p.label = "ERROR"
        return False


def will_be_planet(p, a):
    """
    If the particle's periapsis is within the radius of the planet.
    :param p:
    :param a:
    :return:
    """
    try:
        if p.eccentricity < 1.0 and abs(p.periapsis) <= a:
            p.label = "PLANET"
            return True
        return False
    except Exception as e:
        p.label = "ERROR"
        return False


def is_disk(p, a):
    try:
        if p.eccentricity <= 1.0 and abs(p.periapsis) > a:
            p.label = "DISK"
            return True
        return False
    except Exception as e:
        p.label = "ERROR"
        return False


def is_escape(p, a):
    try:
        if p.eccentricity > 1.0:
            p.label = "ESCAPE"
            return True
        return False
    except Exception as e:
        p.label = "ERROR"
        return False


def log(iteration, error, a,
        NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
        NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING, NEW_MASS_PROTOPLANET, NEW_MASS_DISK, NEW_MASS_ESCAPED,
        total_angular_momentum, planet_density, NUM_PARTICLES_NO_CLASSIFICATION, TOTAL_PARTICLES,
        PARTICLES_BEYOND_ROCHE, MASS_BEYOND_ROCHE, satellite_mass, disk_angular_momentum,
        iron_disk_mass_fraction, iron_disk_mass_fraction_beyond_roche):
    EARTH_MASS = 5.972 * 10 ** 24
    LUNAR_MASS = 7.34767309 * 10 ** 22
    L_EM = 3.5 * 10 ** 34
    print(
        "ITERATION: {}\n"
        "ERROR: {}\n"
        "NEW RADIUS: {} km\n"
        "AVG PLANET DENSITY: {}\n"
        "NUM_PARTICLES_PLANET: {}\n"
        "NUM_PARTICLES_IN_DISK: {}\n"
        "NUM_PARTICLES_ESCAPING: {}\n"
        "NUM_PARTICLES_ERROR: {}\n"
        "TOTAL_PARTICLES: {}\n"
        "NUM_DISK_PARTICLES_BEYOND_ROCHE: {}\n"
        "DISK_MASS_BEYOND_ROCHE: {} M_L (target: 0.92 M_L)\n".format(iteration, error, a / 1000.0, planet_density,
                                                                     NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
                                                                     NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING,
                                                                     NUM_PARTICLES_NO_CLASSIFICATION, TOTAL_PARTICLES,
                                                                     PARTICLES_BEYOND_ROCHE,
                                                                     MASS_BEYOND_ROCHE / LUNAR_MASS)
    )
    print(
        "PROTOPLANET MASS: {} M_E ({} KG) (target: 1 M_E)\n"
        "DISK MASS: {} M_L ({} KG) (target: 1.62 M_L)\n"
        "ESCAPING MASS: {} M_L ({} KG) (target: 0.41 M_L)\n"
        "PREDICTED SATELLITE MASS: {} M_L ({} KG) (target: 1 M_L)".format(
            NEW_MASS_PROTOPLANET / EARTH_MASS, NEW_MASS_PROTOPLANET,
            NEW_MASS_DISK / LUNAR_MASS, NEW_MASS_DISK,
            NEW_MASS_ESCAPED / LUNAR_MASS, NEW_MASS_ESCAPED,
            satellite_mass / LUNAR_MASS, satellite_mass
        )
    )
    print(
        "TOTAL ANGULAR MOMENTUM: {} L_EM ({})\n"
        "DISK ANGULAR MOMEMENTUM: {} L_EM ({}) (target: 0.31 LEM)\n".format(total_angular_momentum / L_EM,
                                                                            total_angular_momentum,
                                                                            disk_angular_momentum / L_EM,
                                                                            disk_angular_momentum)
    )
    print(
        "IRON DISK MASS FRACTION: {}\n"
        "IRON DISK MASS FRACTION BEYOND ROCHE: {}\n\n".format(iron_disk_mass_fraction, iron_disk_mass_fraction_beyond_roche)
    )
