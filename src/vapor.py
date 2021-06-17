import pandas as pd

from src.interpolation import NearestNeighbor1D


def calc_vapor_mass_fraction(particles, phase_path="src/phase_data/duniteS_vapour_curve.txt"):
    phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    nearest_neighbor = NearestNeighbor1D()
    num_particles = 0
    vapor_mass_fraction = 0
    for p in particles:
        if p.label == "DISK":
            num_particles += 1
            entropy_i = p.entropy
            temperature_i = p.temperature
            # entropy_liq = interpolate1d(val=temperature_i, val_array=self.phase_df['temperature'],
            #                             interp_array=self.phase_df['entropy_sol_liq'])
            # entropy_vap = interpolate1d(val=temperature_i, val_array=self.phase_df['temperature'],
            #                             interp_array=self.phase_df['entropy_vap'])
            nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=temperature_i,
                                                                        array=phase_df['temperature'])
            entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
            entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
            if entropy_i < entropy_liq:
                vapor_mass_fraction += 0.0
            elif entropy_liq <= entropy_i <= entropy_vap:
                vapor_mass_fraction += (entropy_i - entropy_liq) / (entropy_vap - entropy_liq)
            elif entropy_i > entropy_vap:
                vapor_mass_fraction += 1.0

    vapor_mass_fraction = vapor_mass_fraction / num_particles

    return vapor_mass_fraction
