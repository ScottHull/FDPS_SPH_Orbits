from math import pi, sqrt, acos, sin, cos, exp
import numpy as np
from copy import copy


class Particle:

    def __init__(self, particle_id, tag, mass, position, velocity, density, internal_energy, pressure,
                 potential_energy, entropy, temperature, mass_grav_body, relative_velocity, calculate_elements=True):
        self.particle_id = particle_id
        self.tag = tag
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.relative_velocity = relative_velocity
        self.density = density
        self.internal_energy = internal_energy
        self.pressure = pressure
        self.potential_energy = potential_energy
        self.entropy = entropy
        self.temperature = temperature
        # self.soundsped = soundspeed

        self.a = None
        self.b = None
        self.avg_planet_density = None

        self.label = None  # PLANET, DISK, or ESCAPE
        self.__G = 6.674 * 10 ** -11
        self.mass_grav_body = mass_grav_body
        self.distance = np.linalg.norm(self.position)  # radial distance from center
        if calculate_elements:
            self.recalculate_elements(mass_grav_body=self.mass_grav_body)

    def __total_momentum_vector(self):
        m_x = self.mass * self.relative_velocity[0]
        m_y = self.mass * self.relative_velocity[1]
        m_z = self.mass * self.relative_velocity[2]
        return m_x, m_y, m_z

    def __angular_momentum_vector(self):
        am = self.mass * np.cross(self.position, self.relative_velocity)
        return am

    def __angular_momentum(self):
        am = np.linalg.norm(self.__angular_momentum_vector())
        return am

    def __node_vector(self):
        return np.cross([0, 0, self.position[0]], self.angular_momentum_vector)

    def __total_orbital_energy(self):
        # kinetic energy, KE = 1/2 m v^2
        self.kinetic_energy = (1.0 / 2.0) * self.mass * (np.linalg.norm(self.relative_velocity) ** 2)
        # vectorized gravitational potential energy, PE = (G M_1 M_2) / r
        self.potential_energy = - (self.__G * self.mass_grav_body * self.mass) / np.linalg.norm(self.position)
        return self.kinetic_energy + self.potential_energy

    def __eccentricity(self):
        try:
            self.alpha = - self.__G * self.mass * self.mass_grav_body
            self.mass_reduced = (self.mass * self.mass_grav_body) / (self.mass + self.mass_grav_body)
            sp_mom = np.cross(self.position, self.relative_velocity)
            ecc = sqrt(1.0 + ((2.0 * self.orbital_energy * (self.angular_momentum ** 2)) / (
                    self.mass_reduced * (self.alpha ** 2))))
            ecc_check = sqrt(1 + 2 * self.__total_orbital_energy() * sp_mom[
                2] ** 2 / self.mass / self.__G / self.__G / self.mass_grav_body / self.mass_grav_body)
            return ecc
        except:
            print("error for particle: {} (ORBITAL ENERGY: {}, ANGULAR MOMENTUM: {})".format(self.particle_id,
                                                                                             self.orbital_energy,
                                                                                             self.angular_momentum))
            return 0

    def __eccentricity_vector(self):
        mu = self.__G * self.mass_grav_body
        term1 = ((((np.linalg.norm(self.relative_velocity)) ** 2) / mu) - (1.0 / self.distance)) * \
                np.array(self.position)
        term2 = ((np.dot(self.position, self.relative_velocity) / mu)) * np.array(
            self.relative_velocity)
        return term1 - term2

    def __semi_major_axis(self):
        E_spec = self.__total_orbital_energy() / self.mass
        mu = self.__G * self.mass_grav_body
        a = - mu / (2.0 * E_spec)
        # a_check = -self.__G * self.mass_grav_body * self.mass / 2 / self.__total_orbital_energy()
        return a

    def __inclination(self):
        return acos(self.angular_momentum_vector[2] / self.angular_momentum) * (180 / pi)

    def __longitude_of_ascending_node(self):
        return np.cross([0, 0, 1], self.angular_momentum_vector)

    def __argument_of_periapsis(self):
        return acos(np.dot(self.periapsis_node_vector, self.eccentricity_vector) /
                    (np.linalg.norm(self.periapsis_node_vector) * np.linalg.norm(self.eccentricity_vector)))

    def __true_anomaly(self):
        """
        http://www.braeunig.us/space/orbmech.htm
        The true anomaly is the position of the orbiting body of the ellipse at a given time.
        The true anomaly equals the mean anomaly for a circular orbit.
        :return:
        """
        return acos(np.dot(self.eccentricity_vector, self.position) /
                    (np.linalg.norm(self.eccentricity_vector) * np.linalg.norm(self.position)))

    def __eccentric_anomaly(self):
        return acos(self.position[0] / self.semi_major_axis)

    def __mean_anomaly(self):
        E = self.__eccentric_anomaly()
        return E - (self.eccentricity * sin(E))

    def __periapsis(self):
        return self.semi_major_axis * (1.0 - self.eccentricity)

    def __additional_heating_from_orbital_circularization(self):
        """
        From Nakajima & Stevenson 2014.  The additional heating delta E_i for each particle due to the circularization
        of the particle's orbit.  This lowers the total energy of the orbit and therefore must heat the particle.
        The resulting semi-major axis would be a_final = a_i (1 - e_i^2) cos^2 (i).
        The delta E_i is the difference between the energy of the circularized vs. eccentric orbit.
        Note that this is a specific energy because particle mass drops.
        :return:
        """
        grav_coeff = (self.__G * self.mass_grav_body) / (2 * self.semi_major_axis)
        eccentric = 1 / ((1 - self.eccentricity ** 2) * (cos(self.inclination * (pi / 180)) ** 2))
        return grav_coeff * (1 - eccentric)

    def __entropy_gain_from_circularization(self):
        """
        The entropy gain from circularization is given by delta S_circ = delta E_circ / delta T.
        Assume that delta U = delta E.  The entropy gain can be written as: delta S = delta U / T.
        Added a minus sign to flip the sign an get an entropy gain for heating.
        :return:
        """
        return - self.circularization_energy_delta / self.temperature


    def recalculate_elements(self, mass_grav_body):
        try:
            self.mass_grav_body = mass_grav_body
            self.distance = np.linalg.norm(self.position)
            self.angular_momentum_vector = self.__angular_momentum_vector()
            self.angular_momentum = self.__angular_momentum()
            self.momentum_vector = self.__total_momentum_vector()
            self.semi_major_axis = self.__semi_major_axis()
            self.orbital_energy = self.__total_orbital_energy()
            self.eccentricity = self.__eccentricity()
            self.eccentricity_vector = self.__eccentricity_vector()
            self.periapsis_node_vector = self.__node_vector()
            self.inclination = self.__inclination()
            self.longitude_of_ascending_node = self.__longitude_of_ascending_node()
            self.argument_of_periapsis = self.__argument_of_periapsis()
            # self.true_anomaly = self.__true_anomaly()
            self.periapsis = self.__periapsis()
            self.circularization_energy_delta = self.__additional_heating_from_orbital_circularization()
            self.circularization_entropy_delta = self.__entropy_gain_from_circularization()
        except Exception as e:
            self.label = "ERROR"
