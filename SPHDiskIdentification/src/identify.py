import numpy as np
import pandas as pd

from src.elements import (total_orbital_energy, angular_momentum, eccentricity, semi_major_axis,
                          equivalent_circular_semi_major_axis, angular_momentum_vector, periapsis, inclination,
                          circularization_entropy_gain)


class ParticleMap:
    """
    Class for mapping particles to their respective location (planet, disk, or escaping particles.
    """

    def __init__(self, particles: pd.DataFrame, mass_planet: float, equatorial_radius: float, center=True,
                 center_on_iron=True):
        self.center_of_mass = None
        self.center_on_iron = center_on_iron
        self.particles = particles  # the dataframe of particles
        self.mass_planet = mass_planet  # the mass of the main planetary body
        self.equatorial_radius = equatorial_radius  # the equatorial radius of the planet
        self.poloidal_radius = equatorial_radius  # the initial poloidal radius of the planet
        self.oblateness = None  # the oblateness of the planet
        self.center = center  # whether to center the particles around the center of mass of the planet
        self.f = None  # the oblateness factor of the planet
        self.bulk_density = None  # the bulk density of the planet
        self.has_converged = False  # whether the iterative solution has converged
        self.iterations = 0  # the number of iterations to converge
        self.error = 1e99  # the error of the iterative solution

    def check_convergence(self, equatorial_radius, target_error=10 ** -8):
        """
        Check if the iterative solution has converged.
        :param target_error: the target error to converge to.
        :param equatorial_radius: the equatorial radius of the planet.
        :return:
        """
        self.error = np.abs(equatorial_radius - self.equatorial_radius) / self.equatorial_radius
        if self.error <= target_error:
            self.has_converged = True
            return True

    def calculate_center_of_mass(self, particles: pd.DataFrame):
        """
        Calculate the center of mass.
        :return:
        """
        total_mass = particles['mass'].sum()
        x = (particles['mass'] * particles['x']).sum() / total_mass
        y = (particles['mass'] * particles['y']).sum() / total_mass
        z = (particles['mass'] * particles['z']).sum() / total_mass
        return x, y, z

    def planetary_center_of_mass(self, target_iron_id=1):
        """
        Calculate the center of mass (COM) of the planetary body.
        This is best done by taking the COM of the target iron, which doesn't get ejected typically.
        Override this function for really, really exotic impacts.
        :return:
        """
        if self.center_on_iron:
            center_on = self.particles[self.particles['tag'] == target_iron_id]
        else:
            center_on = self.particles
        return self.calculate_center_of_mass(center_on)

    def prepare_particles(self):
        """
        Run this to calculate some stuff before identifying particles that is typically not given by SPH outputs.
        :return:
        """
        if self.center:
            self.center_of_mass = self.planetary_center_of_mass()
            self.particles['x'] -= self.center_of_mass[0]
            self.particles['y'] -= self.center_of_mass[1]
            self.particles['z'] -= self.center_of_mass[2]

        self.particles['position'] = np.linalg.norm(self.particles[['x', 'y', 'z']].values, axis=1)
        self.particles['velocity'] = np.linalg.norm(self.particles[['vx', 'vy', 'vz']].values, axis=1)
        # store the position of the particles as a vector in a new column
        self.particles['label'] = None

    def planet_bulk_density(self):
        """
        Calculate the bulk density of the planet.
        :return:
        """
        return self.mass_planet / (4 / 3 * np.pi * self.equatorial_radius ** 3)

    def is_planet(self, particle):
        """
        Determine if a particle is part of the planet.
        Determines if the position of the particle is currently within the equatorial radius of the planet.
        If so, it is by definition part of the planet.
        :param particle:
        :return:
        """
        if particle['position'] <= self.equatorial_radius:
            self.particles.loc[particle.name, 'label'] = 'PLANET'

    def will_be_planet(self, particle):
        """
        Determine if a particle will become part of the planet.
        Determines if the semi-major axis of the circular orbit of the particle is within the equatorial radius
        of the planet.  If so, the particle will eventually accrete to the planet and become part of it.
        :return:
        """
        if particle['circular semi major axis'] <= self.equatorial_radius:
            self.particles.loc[particle.name, 'label'] = 'PLANET'

    def is_disk_or_escaping(self, particle):
        """
        If the periapsis of the particle is greater than the equatorial radius of the planet and it's not
        on a hyperbolic orbit, then it is part of the disk.
        :param particle:
        :return:
        """
        if particle['periapsis'] > self.equatorial_radius and particle['eccentricity'] <= 1:
            self.particles.loc[particle.name, 'label'] = 'DISK'
        elif particle['eccentricity'] > 1:
            self.particles.loc[particle.name, 'label'] = 'ESCAPE'

    def is_planet_disk_or_escaping(self, particle):
        """
        Determine if a particle is part of the planet, disk, or escaping.
        """
        # self.is_planet(particle)
        # self.will_be_planet(particle)
        # self.is_disk_or_escaping(particle)
        is_planet_condition = particle['position'] <= self.equatorial_radius
        will_be_planet_condition = particle['circular semi major axis'] <= self.equatorial_radius
        is_disk_condition = particle['periapsis'] > self.equatorial_radius and particle['eccentricity'] <= 1
        is_escape_condition = particle['eccentricity'] > 1 and particle['position'] > self.equatorial_radius

    def get_label(self):
        """
        Get the label of the particle, either planet, disk, or escaping.
        Use numpy.where to do this faster so that the operation can be vectorized.
        """
        is_planet = (self.particles['position'] <= self.equatorial_radius) | (
                    self.particles['circular semi major axis'] <= self.equatorial_radius)
        is_disk = (self.particles['periapsis'] > self.equatorial_radius) & (self.particles['eccentricity'] <= 1)
        is_escape = (self.particles['eccentricity'] > 1) & (self.particles['position'] > self.equatorial_radius)
        return np.where(is_planet, 'PLANET', np.where(is_disk, 'DISK', np.where(is_escape, 'ESCAPE', None)))

    def roche_radius(self):
        """
        Get the Roche radius of the planet.
        :return: 
        """
        return 2.9 * self.equatorial_radius

    def calculate_elements(self):
        """
        Calculate all necessary orbital elements in order.
        :return:
        """
        print("Calculating orbital elements...")
        # self.particles['angular momentum vector'] = angular_momentum_vector(self.particles)
        # self.particles['z angular momentum vector'] = z_angular_momentum_vector(self.particles)
        self.particles["Lx"], self.particles["Ly"], self.particles["Lz"] = angular_momentum_vector(self.particles)
        self.particles['angular momentum'] = angular_momentum(self.particles)
        self.particles['orbital energy'] = total_orbital_energy(self.particles, self.mass_planet)
        self.particles['eccentricity'] = eccentricity(self.particles, self.mass_planet)
        self.particles['semi major axis'] = semi_major_axis(self.particles, self.mass_planet)
        self.particles['periapsis'] = periapsis(self.particles)
        self.particles['inclination'] = inclination(self.particles)
        self.particles['circular semi major axis'] = equivalent_circular_semi_major_axis(self.particles,
                                                                                         self.mass_planet)
        self.particles['circularization entropy gain'] = circularization_entropy_gain(self.particles, self.mass_planet)
        print("Calculating orbital elements complete.")

    def calculate_planetary_oblateness(self, K=0.335):
        """
        Calculate the oblateness of the planet.
        :param K:
        :return:
        """
        print("Calculating planetary oblateness...")
        G = 6.67408e-11  # m^3 kg^-1 s^-2, gravitational constant
        # sum the third element of the angular momentum vector for all particles
        z_angular_momentum_planet = self.particles['Lz'].sum()
        moment_of_inertia_planet = (2.0 / 5.0) * self.mass_planet * self.equatorial_radius ** 2
        angular_velocity_protoplanet = z_angular_momentum_planet / moment_of_inertia_planet
        keplerian_velocity_protoplanet = np.sqrt(G * self.mass_planet / self.equatorial_radius ** 3)
        numerator = (5.0 / 2.0) * ((angular_velocity_protoplanet / keplerian_velocity_protoplanet) ** 2)
        denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
        f = numerator / denominator
        print(f"Calculating planetary oblateness complete ({f}).")
        return f

    def calculate_planetary_radii(self, mass_planet):
        """
        Calculate the equatorial and polar radii of the planet.
        :param mass_planet:
        :return:
        """
        print("Calculating planetary radii...")
        self.oblateness = self.calculate_planetary_oblateness()
        r = ((3 * mass_planet) / (4 * np.pi * self.bulk_density * (1 - self.oblateness))) ** (1 / 3)
        print(f"Calculating planetary radii complete ({r / 1000} km).")
        return r

    def calculate_planet_mass(self):
        """
        Calculate the mass of the planet.
        :return:
        """
        print("Calculating planetary mass...")
        m = self.particles[self.particles['label'] == 'PLANET']['mass'].sum()
        print(f"Calculating planetary mass complete ({m} kg).")
        return m

    def get_number_of_particles_in_planet_disk_escape(self):
        """
        Get the number of particles in the planet, disk, and escaping.
        :return:
        """
        print("Getting number of particles in planet, disk, and escaping...")
        planet = len(self.particles[self.particles['label'] == 'PLANET'])
        disk = len(self.particles[self.particles['label'] == 'DISK'])
        escape = len(self.particles[self.particles['label'] == 'ESCAPE'])
        not_labelled = len(self.particles[self.particles['label'] == None])
        print(f"Getting number of particles in planet, disk, and escaping complete ({planet}, {disk}, {escape} "
              f"(errors, {not_labelled})).")
        return planet, disk, escape

    def identify(self):
        """
        The main function for identifying particles.
        :return:
        """
        while not self.has_converged:
            print(f"Beginning convergence iteration {self.iterations}")
            # calculate orbital elements
            self.calculate_elements()
            # map the particles to their respective location
            print("Mapping particles to their respective location...")
            # self.particles['label'].apply(self.is_planet_disk_or_escaping, axis=1)
            self.particles['label'] = self.get_label()
            print("Mapping particles to their respective location complete.")
            # calculate the new oblateness, planet mass, and equatorial radius
            mass_planet = self.calculate_planet_mass()
            equatorial_radius = self.calculate_planetary_radii(mass_planet)
            self.error = np.abs(equatorial_radius - self.equatorial_radius) / self.equatorial_radius
            # check for solution convergence
            self.has_converged = self.check_convergence(equatorial_radius)
            # update the planet mass and equatorial radius
            self.mass_planet = mass_planet
            self.equatorial_radius = equatorial_radius
            self.poloidal_radius = self.poloidal_radius * (1 - self.oblateness)
            # get the number of particles in the planet, disk, and escaping
            planet, disk, escape = self.get_number_of_particles_in_planet_disk_escape()
            print(f"Convergence iteration {self.iterations} complete (error: {self.error}).")
            self.iterations += 1

    def loop(self, num_iterations=2):
        """
        Loop over the identification process until convergence.
        :param num_iterations:
        :return:
        """
        self.prepare_particles()  # prepare the particles for identification
        for i in range(num_iterations):
            self.bulk_density = self.planet_bulk_density()
            self.identify()
        return self.particles
