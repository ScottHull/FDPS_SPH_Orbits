import os
import pandas as pd
from math import pi, sqrt
import numpy as np
from copy import copy
import csv
from statistics import mean

from src import elements, centering, classify, vapor


class ParticleMap:

    def __init__(self, path, center, relative_velocity=False, centering_resolution=1e5,
                 centering_delta=1e7, a=(12713.6 / 2.0) * 1000.0, b=(12756.2 / 2.0) * 1000.0, formatted=False):
        self.path = path
        self.time = None
        self.formatted = formatted
        with open(path, 'r') as infile:  # get simulation time from combined file
            self.time = float(infile.readline().strip())
            infile.close()
        if not formatted:
            self.output = pd.read_csv(self.path, skiprows=2, header=None, delimiter="\t")
        else:
            self.output = pd.read_csv(self.path, skiprows=2, delimiter=",")
        # self.a = a  # present-day equatorial radius of the Earth in m
        self.a = 1.099e7
        self.b = b  # present-day polar radius of the Earth in m
        # self.mass_protoearth = classify.calc_mass_protoearth(a=self.a, b=self.b)
        self.mass_protoearth = 5.972e24
        self.center = center
        self.centering_resolution = centering_resolution
        self.centering_delta = centering_delta
        self.com = [0, 0, 0]
        self.vmf = None
        if center:
            if not formatted:
                self.com = centering.center_of_mass(
                    x_coords=self.output[3],
                    y_coords=self.output[4],
                    z_coords=self.output[5],
                    masses=self.output[2],
                    particle_tags=self.output[1],
                    target_iron=True
                )  # COM of the target iron
            else:
                self.com = centering.center_of_mass(
                    x_coords=self.output['x'],
                    y_coords=self.output['y'],
                    z_coords=self.output['z'],
                    masses=self.output['mass'],
                    particle_tags=self.output['tag'],
                    target_iron=True
                )  # COM of the target iron
        self.relative_velocity = relative_velocity
        self.EARTH_MASS = 5.972 * 10 ** 24
        self.LUNAR_MASS = 7.34767309 * 10 ** 22
        self.L_EM = 3.5 * 10 ** 34
        self.targets = {
            "disk mass beyond roche": "0.92 M_L = {} kg".format(0.92 / self.LUNAR_MASS),
            "protoplanet mass": "1 M_E = {} kg".format(1 * self.EARTH_MASS),
            "disk mass": "1.62 M_L = {} kg".format(1.62 * self.LUNAR_MASS),
            "escaping mass": "0.41 M_L = {} kg".format(0.41 * self.LUNAR_MASS),
            "disk angular momentum": "0.31 L_EM = {} kg m^2 / s".format(0.31 * self.L_EM),
            "total angular momentum": "1 L_EM = {} kg m^2 / s".format(1 * self.L_EM),
        }

    def report(self, name="disk_profile.txt"):
        if name in os.listdir(os.getcwd()):
            os.remove(name)
        outfile = open(name, "w")
        for key in self.profile_report.keys():
            value = self.profile_report[key]
            line = "{}\t{}\n".format(key, value)
            outfile.write(line)
        outfile.close()

    def collect_particles(self, find_orbital_elements=True):
        self.particles = classify.collect_particles(
            output=self.output,
            com=self.com,
            mass_protoearth=self.mass_protoearth,
            find_orbital_elements=find_orbital_elements,
            relative_velocity=self.relative_velocity,
            formatted=self.formatted
        )
        return self.particles

    def __convergence_loop(self, particles, K, G):
        iteration = 0
        CONVERGENCE = False
        self.avg_density = classify.average_density(planet_mass=self.mass_protoearth, a=self.a)
        # print("INITIAL AVG DENSITY: {} (iteration: {})".format(self.avg_density, iteration))
        while CONVERGENCE is False:
            NUM_PARTICLES_PLANET = 0
            NUM_PARTICLES_IN_DISK = 0
            NUM_PARTICLES_ESCAPING = 0
            NUM_PARTICLES_NO_CLASSIFICATION = 0
            NEW_MASS_PROTOPLANET = 0.0
            NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET = 0.0
            NEW_MASS_DISK = 0.0
            NEW_Z_ANGULAR_MOMENTUM_DISK = 0.0
            NEW_MASS_ESCAPED = 0.0
            NEW_Z_ANGULAR_MOMENTUM_ESCAPED = 0.0
            PARTICLES_BEYOND_ROCHE = 0
            MASS_BEYOND_ROCHE = 0.0
            for p in particles:
                if classify.is_planet(p=p, a=self.a) or classify.will_be_planet(p=p, a=self.a):
                    NUM_PARTICLES_PLANET += 1
                    NEW_MASS_PROTOPLANET += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                elif classify.is_disk(p=p, a=self.a):
                    NUM_PARTICLES_IN_DISK += 1
                    NEW_MASS_DISK += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_DISK += p.angular_momentum_vector[
                        2]  # assume z component dominate and x and y cancel
                    if classify.is_beyond_roche_radius(p=p):
                        PARTICLES_BEYOND_ROCHE += 1
                        MASS_BEYOND_ROCHE += p.mass
                elif classify.is_escape(p=p, a=self.a):
                    NUM_PARTICLES_ESCAPING += 1
                    NEW_MASS_ESCAPED += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_ESCAPED += p.angular_momentum_vector[2]
                else:
                    NUM_PARTICLES_NO_CLASSIFICATION += 1

            TOTAL_PARTICLES = NUM_PARTICLES_PLANET + NUM_PARTICLES_IN_DISK + NUM_PARTICLES_ESCAPING + \
                              NUM_PARTICLES_NO_CLASSIFICATION
            # recalibrate the system
            iron_disk_mass_fraction, iron_disk_mass_fraction_beyond_roche = classify.get_iron_fraction(
                particles=particles)
            moment_of_inertia_protoplanet = (2.0 / 5.0) * NEW_MASS_PROTOPLANET * (self.a ** 2)
            angular_velocity_protoplanet = NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET / moment_of_inertia_protoplanet
            keplerian_velocity_protoplanet = sqrt((G * NEW_MASS_PROTOPLANET) / (self.a ** 3))
            f_numerator = (5.0 / 2.0) * ((angular_velocity_protoplanet / keplerian_velocity_protoplanet) ** 2)
            f_denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
            new_f = f_numerator / f_denominator
            new_a = ((3 * NEW_MASS_PROTOPLANET) / (4 * pi * self.avg_density * (1 - new_f))) ** (1 / 3)
            error = abs((new_a - self.a) / self.a)
            if error < 10 ** -8:
                CONVERGENCE = True
            else:
                CONVERGENCE = False
            self.a = copy(new_a)
            self.b = self.a * (1 - new_f)
            self.mass_protoearth = copy(NEW_MASS_PROTOPLANET)
            iteration += 1
            total_angular_momentum = sum([i.angular_momentum_vector[2] for i in particles])  # in z-direction
            disk_angular_momentum = sum([np.linalg.norm(i.angular_momentum_vector) for i in particles if i.label == "DISK"])
            satellite_mass = classify.predicted_satellite_mass(
                disk_angular_momentum=NEW_Z_ANGULAR_MOMENTUM_DISK,
                mass_target=NEW_MASS_PROTOPLANET,
                mass_disk=NEW_MASS_DISK,
                mass_escape=NEW_MASS_ESCAPED
            )
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
                    p.a = self.a
                    p.b = self.b
                    p.avg_planet_density = self.avg_density
                    p.recalculate_elements(mass_grav_body=self.mass_protoearth)
                except:
                    pass
            classify.log(
                iteration, error, self.a,
                NUM_PARTICLES_PLANET,
                NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING, NEW_MASS_PROTOPLANET, NEW_MASS_DISK, NEW_MASS_ESCAPED,
                total_angular_momentum, self.avg_density, NUM_PARTICLES_NO_CLASSIFICATION, TOTAL_PARTICLES,
                PARTICLES_BEYOND_ROCHE, MASS_BEYOND_ROCHE, satellite_mass, NEW_Z_ANGULAR_MOMENTUM_DISK,
                iron_disk_mass_fraction, iron_disk_mass_fraction_beyond_roche
            )
            self.profile_report = {
                "num_particles_planet": NUM_PARTICLES_PLANET,
                "num_particles_disk": NUM_PARTICLES_IN_DISK,
                "num_particles_escaping": NUM_PARTICLES_ESCAPING,
                "num_particles_error": NUM_PARTICLES_NO_CLASSIFICATION,
                "total_particles": TOTAL_PARTICLES,
                "particles_beyond_roche": PARTICLES_BEYOND_ROCHE,
                "mass_protoplanet": NEW_MASS_PROTOPLANET,
                "mass_protoplanet (M_E)": NEW_MASS_PROTOPLANET / self.EARTH_MASS,
                "mass_disk": NEW_MASS_DISK,
                "mass_disk (M_L)": NEW_MASS_DISK / self.LUNAR_MASS,
                "mass_escaped": NEW_MASS_ESCAPED,
                "mass_escaped (M_L)": NEW_MASS_ESCAPED / self.LUNAR_MASS,
                "total_angular_momentum": total_angular_momentum,
                "total_angular_momentum (L_EM)": total_angular_momentum / self.L_EM,
                "disk_angular_momentum": disk_angular_momentum,
                "disk_angular_momentum (L_EM)": disk_angular_momentum / self.L_EM,
                "z_angular_momentum_disk": NEW_Z_ANGULAR_MOMENTUM_DISK,
                "z_angular_momentum_disk (L_EM)": NEW_Z_ANGULAR_MOMENTUM_DISK / self.L_EM,
                "average_density": self.avg_density,
                "disk_mass_beyond_roche": MASS_BEYOND_ROCHE,
                "disk_mass_beyond_roche (M_L)": MASS_BEYOND_ROCHE / self.LUNAR_MASS,
                "satellite_mass": satellite_mass,
                "iron_disk_mass_fraction": iron_disk_mass_fraction,
                "iron_disk_mass_fraction_beyond_roche": iron_disk_mass_fraction_beyond_roche,
                "a": self.a,
                "b": self.b,
            }
            self.num_particles_planet = NUM_PARTICLES_PLANET
            self.num_particles_disk = NUM_PARTICLES_IN_DISK
            self.num_particles_escaping = NUM_PARTICLES_ESCAPING

    def solve(self, particles, phase_path, K=0.335, G=6.674 * 10 ** -11, report_name="disk_profile.txt", iteration="UNSPECIFIED", simulation_time="UNSPECIFIED"):
        # K = 0.335 for Earth, K = 2/5 for homogenous body
        self.__convergence_loop(particles=particles, K=K, G=G)
        self.__convergence_loop(particles=particles, K=K, G=G)  # run twice to recalc avg density after initial solution
        if phase_path is None:
            self.vmf = None
        else:
            self.vmf = vapor.calc_vapor_mass_fraction(particles=particles, phase_path=phase_path)
        self.profile_report.update({"disk vmf": self.vmf})
        self.profile_report.update({"name": report_name})
        self.profile_report.update({'phase_path': phase_path})
        self.profile_report.update({'iteration profiled (corresponds to output number)': iteration})
        self.profile_report.update({'simulation time': simulation_time})
        self.report(name=report_name)
