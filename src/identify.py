import pandas as pd
from math import pi, sqrt

from src import elements, centering, classify


class ParticleMap:

    def __init__(self, path, center, relative_velocity=False, centering_resolution=1e5,
                 centering_delta=1e7, a=(12713.6 / 2.0) * 1000.0, b=(12756.2 / 2.0) * 1000.0):
        self.path = path
        self.output = pd.read_csv(self.path, skiprows=2, header=None, delimiter="\t")
        self.a = a  # present-day equatorial radius of the Earth in m
        self.b = b  # present-day polar radius of the Earth in m
        self.mass_protoearth = classify.calc_mass_protoearth(a=self.a, b=self.b)
        self.center = center
        self.centering_resolution = centering_resolution
        self.centering_delta = centering_delta
        self.com = [0, 0, 0]
        if center:
            self.com = centering.center_of_mass(
                x_coords=self.output[3],
                y_coords=self.output[4],
                z_coords=self.output[5],
                masses=self.output[2],
                particle_tags=self.output[1],
                target_iron=True
            )  # COM of the target iron
        self.relative_velocity = relative_velocity

    def collect_particles(self, find_orbital_elements=True):
        return classify.collect_particles(
            output=self.output,
            com=self.com,
            mass_protoearth=self.mass_protoearth,
            find_orbital_elements=find_orbital_elements
        )

    def solve(self, particles, K=0.335, G=6.674 * 10 ** -11, avg_density=5.5 * 1000):
        iteration = 0
        CONVERGENCE = False
        while CONVERGENCE is False:
            NUM_PARTICLES_WITHIN_RADIAL_DISTANCE = 0
            NUM_PARTICLES_IN_DISK = 0
            NUM_PARTICLES_ESCAPING = 0
            NUM_PARTICLES_NO_CLASSIFICATION = 0
            NEW_MASS_PROTOPLANET = 0.0
            NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET = 0.0
            NEW_MASS_DISK = 0.0
            NEW_Z_ANGULAR_MOMENTUM_DISK = 0.0
            NEW_MASS_ESCAPED = 0.0
            NEW_Z_ANGULAR_MOMENTUM_ESCAPED = 0.0
            for p in particles:
                if classify.is_planet(p=p, a=self.a) or classify.will_be_planet(p=p, a=self.a):
                    NUM_PARTICLES_WITHIN_RADIAL_DISTANCE += 1
                    NEW_MASS_PROTOPLANET += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                elif classify.is_disk(p=p, a=self.a):
                    NUM_PARTICLES_IN_DISK += 1
                    NEW_MASS_DISK += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_DISK += p.angular_momentum_vector[
                        2]  # assume z component dominate and x and y cancel
                elif classify.is_escape(p=p, a=self.a):
                    NUM_PARTICLES_ESCAPING += 1
                    NEW_MASS_ESCAPED += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_ESCAPED += p.angular_momentum_vector[2]
                else:
                    NUM_PARTICLES_NO_CLASSIFICATION += 1

            # recalibrate the system
            moment_of_inertia_protoplanet = (2.0 / 5.0) * NEW_MASS_PROTOPLANET * (self.a ** 2)
            angular_velocity_protoplanet = NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET / moment_of_inertia_protoplanet
            keplerian_velocity_protoplanet = sqrt((G * NEW_MASS_PROTOPLANET) / self.a ** 3)
            f_numerator = (5.0 / 2.0) * ((angular_velocity_protoplanet / keplerian_velocity_protoplanet) ** 2)
            f_denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
            new_f = f_numerator / f_denominator
            new_a = ((3.0 * pi * NEW_MASS_PROTOPLANET * (1.0 - new_f)) / (4.0 * avg_density)) ** (1 / 3)
            error = abs((new_a - self.a) / self.a)
            if error < 10 ** -8:
                CONVERGENCE = True
            else:
                CONVERGENCE = False
            self.a = new_a
            self.mass_protoearth = NEW_MASS_PROTOPLANET
            iteration += 1
            total_angular_momentum = sum([i.angular_momentum for i in particles])
            if self.relative_velocity:
                new_target_velocity = classify.refine_target_velocity(particles=particles)
                for p in particles:
                    if self.relative_velocity is True:
                        p.relative_velocity_vector = [
                            p.velocity[0] - new_target_velocity[0],
                            p.velocity[1] - new_target_velocity[1],
                            p.velocity[2] - new_target_velocity[2]
                        ]
            for p in particles:
                try:
                    p.recalculate_elements(mass_grav_body=self.mass_protoearth)
                except:
                    pass
            classify.log(
                iteration, error, self.a,
                NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
                NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING, NEW_MASS_PROTOPLANET, NEW_MASS_DISK, NEW_MASS_ESCAPED,
                total_angular_momentum
            )


class ParticleMapFromFiles:

    def __init__(self, path):
        self.path = path

    def read(self, time):
        particles = []
        df = pd.read_csv(self.path + "/{}.csv".format(time))
