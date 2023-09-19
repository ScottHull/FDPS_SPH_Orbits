import numpy as np
import pandas as pd

G = 6.674 * 10 ** -11  # m^3 kg^-1 s^-2, the gravitational constant


def angular_momentum_vector(df: pd.DataFrame):
    """
    Returns a 3D vector of the angular momentum of a particle.
    """
    position = np.array(df[['x', 'y', 'z']].values)
    velocity = np.array(df[['vx', 'vy', 'vz']].values)
    # return the angular momentum vector as lists of the x, y, and z components
    return (df['mass'][:, np.newaxis] * np.cross(position, velocity)).T.tolist()


def z_angular_momentum_vector(df: pd.DataFrame):
    """
    Return the z component of the angular momentum of a particle.
    """
    return df['angular momentum vector'][2]


def angular_momentum(df: pd.DataFrame):
    """
    Returns the magnitude of the angular momentum of a particle.
    """
    # return np.linalg.norm(df['angular momentum vector'].tolist())
    angular_momentum_components = df[['Lx', 'Ly', 'Lz']].values
    return np.linalg.norm(angular_momentum_components, axis=1)


def total_orbital_energy(df: pd.DataFrame, mass_grav_body: float):
    """
    Returns the total orbital energy of a particle.
    """
    # kinetic energy, KE = 1/2 m v^2
    kinetic_energy = 0.5 * df['mass'] * (df['velocity'] ** 2)
    # vectorized gravitational potential energy, PE = (G M_1 M_2) / r
    potential_energy = - (G * mass_grav_body * df['mass']) / df['position']
    return kinetic_energy + potential_energy


def eccentricity(df: pd.DataFrame, mass_grav_body: float):
    """
    Returns the eccentricity of a particle.
    """
    alpha = - G * df['mass'] * mass_grav_body
    reduced_mass = (df['mass'] * mass_grav_body) / (df['mass'] + mass_grav_body)
    return np.sqrt(1 + ((2 * df['orbital energy'] * (df['angular momentum'] ** 2)) / (reduced_mass * alpha ** 2)))


def semi_major_axis(df: pd.DataFrame, mass_grav_body: float):
    """
    Returns the semi-major axis of a particle.
    """
    return - (G * mass_grav_body) / (2 * df['orbital energy'] / df['mass'])


def inclination(df: pd.DataFrame):
    """
    Returns the inclination of a particle.
    """
    return np.arccos(df['Lz'] / df['angular momentum']) * (180 / np.pi)


def periapsis(df: pd.DataFrame):
    """
    Returns the periapsis of a particle.
    """
    return df['semi major axis'] * (1.0 - df['eccentricity'])


def equivalent_circular_semi_major_axis(df: pd.DataFrame, mass_grav_body: float):
    """
    Many SPH studies use the equivalent circular semi-major axis to determine if a particle is in the disk.
    The equivalent circular orbit is a circularized (i.e., no eccentricity) orbit with the same angular momentum as
    the non-circularized orbit. This function returns the equivalent circular semi-major axis of a particle.
    """
    specific_angular_momentum = df['angular momentum'] / df['mass']
    return specific_angular_momentum ** 2 / (G * mass_grav_body)

def circularization_energy_gain(df: pd.DataFrame, mass_grav_body: float):
    """
    The internal energy gain associated with orbital circularization of a particle.
    """
    grav_coeff = G * mass_grav_body / (2 * df['semi major axis'])
    eccentric = 1 / ((1 - df['eccentricity'] ** 2) * (np.cos(df['inclination'] * (np.pi / 180)) ** 2))
    return grav_coeff * (1 - eccentric)

def circularization_entropy_gain(df: pd.DataFrame, mass_grav_body: float):
    """
    Returns the circularization entropy gain of a particle.
    """
    energy_gain = circularization_energy_gain(df, mass_grav_body)
    return - energy_gain / df['temperature']
