import os
import pandas as pd
import numpy as np
from random import randint
from statistics import mean
import multiprocessing.dummy as mp

from src.combine import CombineFile
from src.identify import ParticleMap
from src.vapor import calc_vapor_mass_fraction_with_circularization, calc_vapor_mass_fraction_from_formatted

MASS_EARTH = 5.972 * 10 ** 24
MASS_MOON = 7.34767309 * 10 ** 22
L_EM = 3.5 * 10 ** 34
ROCHE_LIM = 2.9 * (6371 * 1000)
G = 6.67408 * 10 ** -11


def predicted_moon_mass(disk):
    term1 = 1.9 * (sum(disk['angular_momentum']) / ((G * MASS_EARTH * ROCHE_LIM) ** (1 / 2)))
    term2 = 1.1 * sum(disk['mass'])
    term3 = 1.9 * (0.05 * sum(disk['mass']))  # assumption made by Canup & Asphaug 2001, escpaing mass = 0.05 M_D
    return (term1 - term2 - term3) / MASS_MOON


def __mean(l):
    try:
        return mean(l)
    except:
        return 0.0


def get_sim_report(particle_df, phase_path, to_path, iteration, formatted_time, sim_name):
    planet = particle_df[particle_df['label'] == "PLANET"]
    disk = particle_df[particle_df['label'] == "DISK"]
    filtered_disk = disk[disk['circ_entropy_delta'] <= 5000]  # have some highly inclined orbits sometimes
    escaping = particle_df[particle_df['label'] == "ESCAPE"]
    disk_iron = disk[disk['tag'] % 2 != 0]
    from_target = disk[disk['tag'] <= 1]
    from_theia = disk[disk['tag'] > 1]

    disk_iron_mass_fraction = 0.0
    disk_iron_mass_fraction_beyond_roche = 0.0
    theia_disk_mass_fraction = 0.0
    try:
        disk_iron_mass_fraction = sum(disk_iron['mass']) / sum(disk['mass']) * 100.0
    except:
        pass
    try:
        disk_iron_mass_fraction_beyond_roche = sum(disk_iron[disk_iron['radius'] > ROCHE_LIM]['mass']) / sum(
            disk['mass']) * 100.0
    except:
        pass
    try:
        theia_disk_mass_fraction = (sum(from_theia['mass']) / (
                    sum(from_theia['mass']) + sum(from_target['mass']))) * 100.0
    except:
        pass

    vmf = calc_vapor_mass_fraction_with_circularization(particles=particle_df, phase_path=phase_path)
    vmf_without_circ = calc_vapor_mass_fraction_from_formatted(df=particle_df, phase_path=phase_path)
    data = {
        "NAME": sim_name,
        "ITERATION": iteration,
        "TIME_HRS": formatted_time,
        "PLANET_MASS": "{} M_E".format(sum(planet['mass']) / MASS_EARTH),
        "DISK_MASS": "{} M_L".format(sum(disk['mass']) / MASS_MOON),
        "ESCAPING_MASS": "{} M_L".format(sum(escaping['mass']) / MASS_MOON),
        "NUM_PARTICLES_PLANET": len(planet),
        "NUM_PARTICLES_DISK": len(disk),
        "NUM_PARTICLES_ESCAPING": len(escaping),
        "TOTAL_PARTICLES": len(particle_df),
        "DISK_MASS_BEYOND_ROCHE": "{} M_L".format(sum(disk[disk['radius'] > ROCHE_LIM]['mass']) / MASS_MOON),
        "NUM_PARTICLES_BEYOND_ROCHE": len(disk[disk['radius'] > ROCHE_LIM]),
        "DISK_IRON_MASS_FRACTION": "{} %".format(disk_iron_mass_fraction),
        "DISK_IRON_MASS_FRACTION_BEYOND_ROCHE": "{} %".format(disk_iron_mass_fraction_beyond_roche),
        "AVERAGE_PLANET_DENSITY": "{} kg/m3".format(__mean(planet['planet_avg_density'])),
        "PLANET_EQUATORIAL_RADIUS": "{} km".format(__mean(planet['planet_radius_equatorial']) / 1000),
        "PLANET_POLAR_RADIUS": "{} km".format(__mean(planet['planet_radius_polar']) / 1000),
        "DISK_ANGULAR_MOMENTUM": "{} L_EM".format(sum(disk['angular_momentum']) / L_EM),
        "DISK_ANGULAR_MOMENTUM_BEYOND_ROCHE": "{} L_EM".format(
            sum(disk[disk['radius'] > ROCHE_LIM]['angular_momentum']) / L_EM),
        "DISK_VMF_W_CIRC": "{} %".format(vmf * 100),
        "DISK_VMF_WITHOUT_CIRC": "{} %".format(vmf_without_circ * 100),
        "TOTAL_ANGULAR_MOMENTUM": "{} L_EM".format(sum(particle_df['angular_momentum']) / L_EM),
        "MEAN_DISK_ENTROPY_W_CIRC": __mean(filtered_disk['entropy'] + filtered_disk['circ_entropy_delta']),
        "MEAN_DISK_ENTROPY_WITHOUT_CIRC": __mean(disk['entropy']),
        "DISK_DELTA_S_DUE_TO_ORBIT_CIRCULAR_FILTERED": __mean(filtered_disk['circ_entropy_delta']),
        "PREDICTED_MOON_MASS": "{} M_L".format(predicted_moon_mass(disk)),
        "DISK_THEIA_MASS_FRACTION": "{} %".format(theia_disk_mass_fraction),
        "MEAN_DISK_TEMPERATURE": "{} K".format(__mean(disk['temperature'])),
    }
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    pd.DataFrame(data, index=[0]).to_csv(to_path + "/{}.csv".format(iteration))


def write_report_at_time(particles, fname):
    print("writing report " + fname)
    df = pd.DataFrame({
        "id": [p.particle_id for p in particles],
        "x": [p.position[0] for p in particles],
        "y": [p.position[1] for p in particles],
        "z": [p.position[2] for p in particles],
        "vx": [p.velocity[0] for p in particles],
        "vy": [p.velocity[1] for p in particles],
        "vz": [p.velocity[2] for p in particles],
        "mass": [p.mass for p in particles],
        "radius": [p.distance for p in particles],
        "tag": [p.tag for p in particles],
        "label": [p.label for p in particles],
        "entropy": [p.entropy for p in particles],
        "temperature": [p.temperature for p in particles],
        "density": [p.density for p in particles],
        "pressure": [p.pressure for p in particles],
        "internal_energy": [p.internal_energy for p in particles],
        "potential_energy": [p.potential_energy for p in particles],
        "eccentricity": [p.eccentricity for p in particles],
        "inclination": [p.inclination for p in particles],
        "orbital_energy": [p.orbital_energy for p in particles],
        "semi_major_axis": [p.semi_major_axis for p in particles],
        "mass_grav_body": [p.mass_grav_body for p in particles],
        "circ_energy_delta": [p.circularization_energy_delta for p in particles],
        "circ_entropy_delta": [p.circularization_entropy_delta for p in particles],
        "angular_momentum": [p.angular_momentum for p in particles],
        "planet_radius_equatorial": [p.a for p in particles],
        "planet_radius_polar": [p.b for p in particles],
        "planet_avg_density": [p.avg_planet_density for p in particles],
    })
    df.to_csv(fname)
    return df

rows_map = {
        "PLANET_MASS": "{Planet Mass ($M_\oplus$)}",
        "DISK_MASS": "{Disk Mass ($M_L$)}",
        "ESCAPING_MASS": "{Escaping Mass ($M_L$)}",
        "NUM_PARTICLES_PLANET": "{$N_{planet}$}",
        "NUM_PARTICLES_DISK": "{$N_{disk}$}",
        "NUM_PARTICLES_ESCAPING": "{$N_{escape}$}",
        "TOTAL_PARTICLES": "{$N_{total}$}",
        "DISK_MASS_BEYOND_ROCHE": "{Disk Mass $\geq$ $R_{Roche}$ ($M_{L}$)}",
        "NUM_PARTICLES_BEYOND_ROCHE": "{$N_{disk}$ $\geq$ $R_{Roche}$}",
        "DISK_IRON_MASS_FRACTION": "{Disk Iron Mass Fraction ($\%$)}",
        "DISK_MASS_FRACTION_BEYOND_ROCHE": "{Disk Iron Mass Fraction $\geq$ $R_{Roche}$ ($\%$)}",
        "AVERAGE_PLANET_DENSITY": "{Planet Avg. Density ($kg/m^3$)}",
        "PLANET_EQUATORIAL_RADIUS": "{Planet $a$ (km)}",
        "PLANET_POLAR_RADIUS": "{Planet $b$ (km)}",
        "DISK_ANGULAR_MOMENTUM": "{$L_{Disk}$ ($L_{EM}$)}",
        "DISK_ANGULAR_MOMENTUM_BEYOND_ROCHE": "{$L_{Disk}$ $\geq R_{Roche}$ ($L_{EM}$)}",
        "TOTAL_ANGULAR_MOMENTUM": "{L_{total} ($L_{EM}$)}",
        "DISK_VMF_W_CIRC": "{Disk VMF  ($\%$)}",
        "DISK_VMF_WITHOUT_CIRC": "{Disk VMF  ($\%$)}",
        "MEAN_DISK_ENTROPY_W_CIRC": "{Avg. $S_{disk}$ w/ circ. (J/K)}",
        "MEAN_DISK_ENTROPY_WITHOUT_CIRC": "{Avg. $S_{disk}$ w/o circ. (J/K)}",
        "DISK_DELTA_S_DUE_TO_ORBIT_CIRCULAR_FILTERED": "{Avg. $\Delta S_{circ}$ (J/K)}",
        "PREDICTED_MOON_MASS": "{$M_{M}$ ($M_L$)}",
        "DISK_THEIA_MASS_FRACTION": "{Disk Theia Mass Fraction ($\%$)}",
        "MEAN_DISK_TEMPERATURE": "{Avg. Disk Temperature (K)}"
    }

def build_latex_table_from_disk_report(run_names: list, run_titles: list, to_base_path: str, filename: str, iteration: int):
    table_header = ["Simulation"]
    for index, run in enumerate(run_names):
        table_header.append(run_titles[index])
    run_rows = []
    for header in rows_map.keys():
        row = [rows_map[header]]
        for index, run in enumerate(run_names):
            try:
                path = to_base_path + "{}/{}_reports/".format(run, run)
                df = pd.read_csv(path + "{}.csv".format(iteration))
                val = str(df[header][0]).split(" ")[0]
                if "TOTAL_ANGULAR_MOMENTUM" in header and float(val) > 1e30:
                    val = float(val) / L_EM
                if "particle" in header.lower():
                    row.append(str(int(val)))
                else:
                    row.append(str(round(float(val), 2)))
            except Exception as e:
                print(e)
                row.append("")
        run_rows.append(row)

    if os.path.exists(filename):
        os.remove(filename)

    with open(filename, 'w') as outfile:
        header_line = ""
        for index, h in enumerate(table_header):
            header_line += "{" + h + "}"
            if index != len(table_header) - 1:
                header_line += " & "
        header_line += " \\\ \midrule\n"
        outfile.write(header_line)
        for row in run_rows:
            line = ""
            for index, v in enumerate(row):
                if index == 0:
                    line += v
                else:
                    line += "{" + v + "}"
                if index != len(row) - 1:
                    line += " & "
            line += " \\\ \midrule\n"
            outfile.write(line)


class BuildReports:

    def __init__(self, to_dir, from_dir, start_time, end_time, number_processes, eos_phase_path, accessory_path,
                 interval=1, force_solve_at_time=False):
        self.to_dir = to_dir
        self.from_dir = from_dir
        self.start_time = start_time
        self.end_time = end_time
        self.number_processes = number_processes
        self.eos_phase_path = eos_phase_path
        self.interval = interval
        self.labels = {}
        self.accessory_path = accessory_path
        self.force_solve_at_time = force_solve_at_time
        try:
            os.mkdir(accessory_path)
        except:
            pass
        self.__get_end_state()

    def __output_disk_state(self, particles, time, vmf):
        if self.accessory_path is not None:
            output = open("{}/{}.txt".format(self.accessory_path, time), 'w')
            avg_disk_entropy = [p.entropy for p in particles if p.label == "DISK"]
            try:
                avg_disk_entropy = sum(avg_disk_entropy) / len(avg_disk_entropy)
            except:
                avg_disk_entropy = 0
            disk_mass = sum([p.mass for p in particles if p.label == "DISK"])
            line = "time\t{}\ndisk vmf\t{}\navg disk entropy\t{}\ndisk mass\t{}\n".format(time, vmf, avg_disk_entropy,
                                                                                          disk_mass)
            output.write(line)
            output.close()

    def __build_df_from_endstate(self, particles):
        try:
            end_labels = [self.labels[p.particle_id] for p in particles]
        except:
            end_labels = []
        return pd.DataFrame({
            "id": [p.particle_id for p in particles],
            "x": [p.position[0] for p in particles],
            "y": [p.position[1] for p in particles],
            "z": [p.position[2] for p in particles],
            "vx": [p.velocity[0] for p in particles],
            "vy": [p.velocity[1] for p in particles],
            "vz": [p.velocity[2] for p in particles],
            "mass": [p.mass for p in particles],
            "radius": [p.distance for p in particles],
            "tag": [p.tag for p in particles],
            "end_label": end_labels,
            "label": [p.label for p in particles],
            "entropy": [p.entropy for p in particles],
            "temperature": [p.temperature for p in particles],
            "density": [p.density for p in particles],
            "pressure": [p.pressure for p in particles],
            "internal_energy": [p.internal_energy for p in particles],
            "potential_energy": [p.potential_energy for p in particles],
            "eccentricity": [p.eccentricity for p in particles],
            "inclination": [p.inclination for p in particles],
            "orbital_energy": [p.orbital_energy for p in particles],
            "semi_major_axis": [p.semi_major_axis for p in particles],
            "mass_grav_body": [int(p.mass_grav_body) for p in particles],
        })

    def __get_end_state(self):
        particles = self.build_report_at_time(time=self.end_time, solve=True, save=False)
        self.labels = dict([(p.particle_id, p.label) for p in particles])
        self.build_report_at_time(time=self.end_time, solve=False, save=True, end_state=True)

    def build_report_at_time(self, time, pm=None, formatted_time=None, total_particles=None, solve=False, save=True,
                             end_state=False):
        if self.force_solve_at_time and not end_state:
            solve = True
        to_fname = "{}.csv".format(time)
        if not pm:
            fname = "merged_{}_{}_reportsproc.dat".format(time, randint(0, 1000000))
            cf = CombineFile(num_processes=self.number_processes, time=time, output_path=self.from_dir, to_fname=fname)
            combined_file = cf.combine()
            formatted_time = cf.sim_time
            total_particles = cf.total_particles
            f = os.getcwd() + "/{}".format(fname)
            pm = ParticleMap(path=f, center=True, relative_velocity=False)
            particles = pm.collect_particles(find_orbital_elements=solve)
            if solve:
                pm.solve(particles=particles, phase_path=self.eos_phase_path)
                self.__output_disk_state(time=time, particles=particles, vmf=pm.vmf)
            os.remove(f)
        else:
            particles = pm.particles
        if save:
            particle_df = self.__build_df_from_endstate(particles=particles)
            particle_df.to_csv(self.to_dir + "/{}".format(to_fname))
            with open(self.to_dir + "/{}".format(to_fname), 'r+') as infile:
                content = infile.read()
                infile.seek(0, 0)
                infile.write("{}\n{}\n".format(formatted_time, total_particles) + content)
            infile.close()
        return particles

    def make_reports(self, mp_pool_size=5):
        pool = mp.Pool(mp_pool_size)
        for time in np.arange(self.start_time, self.end_time + self.interval, self.interval):
            pool.map(self.build_report_at_time, [time])
        pool.close()
        pool.join()
