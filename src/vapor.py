import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

from src.interpolation import NearestNeighbor1D


def calc_vapor_mass_fraction(particles, phase_path="src/phase_data/duniteS_vapour_curve.txt"):
    phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    nearest_neighbor = NearestNeighbor1D()
    num_particles = 0
    vapor_mass_fraction = 0
    for p in particles:
        if p.label == "DISK" and p.tag % 2 == 0:
            num_particles += 1
            entropy_i = p.entropy
            temperature_i = p.temperature
            # entropy_liq = interpolate1d(val=temperature_i, val_array=self.phase_df['temperature'],
            #                             interp_array=self.phase_df['entropy_sol_liq'])
            # entropy_vap = interpolate1d(val=temperature_i, val_array=self.phase_df['temperature'],
            #                             interp_array=self.phase_df['entropy_vap'])
            nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=temperature_i,
                                                                        array=list(phase_df['temperature']))
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

def plot_disk_entropy(particles, phase_path="src/phase_data/duniteS_vapour_curve.txt"):
    phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.distance / 1000 for p in particles if p.label == "DISK" and p.tag == 0],
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 0],
        color='red',
        marker="+",
        label="DISK // Target // Silicate",
        alpha=0.6
    )
    # ax.scatter(
    #     [p.distance / 1000 for p in particles if p.label == "DISK" and p.tag == 1],
    #     [p.entropy for p in particles if p.label == "DISK" and p.tag == 1],
    #     color='blue',
    #     marker="+",
    #     label="DISK // Target // Iron",
    #     alpha=0.6
    # )
    ax.scatter(
        [p.distance / 1000 for p in particles if p.label == "DISK" and p.tag == 2],
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 2],
        color='orange',
        marker="*",
        label="DISK // Impactor // Silicate",
        alpha=0.6
    )
    # ax.scatter(
    #     [p.distance / 1000 for p in particles if p.label == "DISK" and p.tag == 3],
    #     [p.entropy for p in particles if p.label == "DISK" and p.tag == 3],
    #     color='pink',
    #     marker="*",
    #     label="DISK // Impactor // Iron",
    #     alpha=0.6
    # )
    ax.set_xlabel("Radial Distance (km)")
    ax.set_ylabel("Entropy")
    ax.set_title("Disk Entropy")
    ax.grid()
    ax.legend()
    plt.savefig("entropy.png", format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.plot(
        phase_df['entropy_sol_liq'],
        phase_df['temperature'],
        label='liquid phase boundary',
        color='blue',
        linewidth=2.0,
    )
    ax.plot(
        phase_df['entropy_vap'],
        phase_df['temperature'],
        label='vapor phase boundary',
        color='red',
        linewidth=2.0
    )
    ax.scatter(
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 0],
        [p.temperature for p in particles if p.label == "DISK" and p.tag == 0],
        color='red',
        marker="+",
        label="DISK // Target // Silicate",
        alpha=0.6
    )
    # ax.scatter(
    #     [p.entropy for p in particles if p.label == "DISK" and p.tag == 1],
    #     [p.temperature for p in particles if p.label == "DISK" and p.tag == 1],
    #     color='blue',
    #     marker="+",
    #     label="DISK // Target // Iron",
    #     alpha=0.6
    # )
    ax.scatter(
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 2],
        [p.temperature for p in particles if p.label == "DISK" and p.tag == 2],
        color='orange',
        marker="*",
        label="DISK // Impactor // Silicate",
        alpha=0.6
    )
    # ax.scatter(
    #     [p.entropy for p in particles if p.label == "DISK" and p.tag == 3],
    #     [p.temperature for p in particles if p.label == "DISK" and p.tag == 3],
    #     color='pink',
    #     marker="*",
    #     label="DISK // Impactor // Iron",
    #     alpha=0.6
    # )
    ax.set_xlabel("Entropy")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Disk Entropy vs. Temperature")
    ax.grid()
    ax.legend()
    plt.savefig("temperature_vs_entropy.png", format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    cmap = plt.get_cmap("viridis")
    scales = [p.density for p in particles if p.label == "DISK"]
    norm = plt.Normalize(min(scales), max(scales))
    ax.plot(
        phase_df['entropy_sol_liq'],
        phase_df['temperature'],
        label='liquid phase boundary',
        color='orange',
        linewidth=2.0,
        alpha=0.8
    )
    ax.plot(
        phase_df['entropy_vap'],
        phase_df['temperature'],
        label='vapor phase boundary',
        color='pink',
        linewidth=2.0,
        alpha=0.8
    )
    ax.scatter(
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 0],
        [p.temperature for p in particles if p.label == "DISK" and p.tag == 0],
        c=[p.density for p in particles if p.label == "DISK" and p.tag == 0],
        marker="+",
        label="DISK // Target // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    ax.scatter(
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 2],
        [p.temperature for p in particles if p.label == "DISK" and p.tag == 2],
        c=[p.density for p in particles if p.label == "DISK" and p.tag == 2],
        marker="*",
        label="DISK // Impactor // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.ax.set_title("Density")
    ax.set_xlabel("Entropy")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Disk Entropy vs. Temperature (Density Colored)")
    ax.grid()
    ax.legend()
    plt.savefig("temperature_vs_entropy_and_density.png", format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    cmap = plt.get_cmap("viridis")
    scales = [p.internal_energy for p in particles if p.label == "DISK"]
    norm = plt.Normalize(min(scales), max(scales))
    ax.scatter(
        [p.density for p in particles if p.label == "DISK" and p.tag == 0],
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 0],
        c=[p.internal_energy for p in particles if p.label == "DISK" and p.tag == 0],
        marker="+",
        label="DISK // Target // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    ax.scatter(
        [p.density for p in particles if p.label == "DISK" and p.tag == 2],
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 2],
        c=[p.internal_energy for p in particles if p.label == "DISK" and p.tag == 2],
        marker="*",
        label="DISK // Impactor // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.ax.set_title("Internal Energy")
    ax.set_xlabel("Density")
    ax.set_ylabel("Entropy")
    ax.set_title("Density vs Entropy (Internal Energy Colored)")
    ax.grid()
    ax.legend()
    plt.savefig("density_vs_entropy.png", format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    cmap = plt.get_cmap("viridis")
    scales = [p.entropy for p in particles if p.label == "DISK"]
    norm = plt.Normalize(min(scales), max(scales))
    ax.scatter(
        [p.distance for p in particles if p.label == "DISK" and p.tag == 0],
        [p.density for p in particles if p.label == "DISK" and p.tag == 0],
        c=[p.entropy for p in particles if p.label == "DISK" and p.tag == 0],
        marker="+",
        label="DISK // Target // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    ax.scatter(
        [p.distance for p in particles if p.label == "DISK" and p.tag == 2],
        [p.density for p in particles if p.label == "DISK" and p.tag == 2],
        c=[p.entropy for p in particles if p.label == "DISK" and p.tag == 2],
        marker="*",
        label="DISK // Impactor // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.ax.set_title("Entropy")
    ax.set_xlabel("Radius")
    ax.set_ylabel("Density")
    ax.set_title("Radius vs Density (Entropy Colored)")
    ax.grid()
    ax.legend()
    plt.savefig("radius_vs_density.png", format='png')

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    cmap = plt.get_cmap("viridis")
    scales = [p.distance for p in particles if p.label == "DISK"]
    norm = plt.Normalize(min(scales), max(scales))
    ax.plot(
        phase_df['entropy_sol_liq'],
        phase_df['temperature'],
        label='liquid phase boundary',
        color='orange',
        linewidth=2.0,
        alpha=0.8
    )
    ax.plot(
        phase_df['entropy_vap'],
        phase_df['temperature'],
        label='vapor phase boundary',
        color='pink',
        linewidth=2.0,
        alpha=0.8
    )
    ax.scatter(
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 0],
        [p.temperature for p in particles if p.label == "DISK" and p.tag == 0],
        c=[p.distance for p in particles if p.label == "DISK" and p.tag == 0],
        marker="+",
        label="DISK // Target // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    ax.scatter(
        [p.entropy for p in particles if p.label == "DISK" and p.tag == 2],
        [p.temperature for p in particles if p.label == "DISK" and p.tag == 2],
        c=[p.distance for p in particles if p.label == "DISK" and p.tag == 2],
        marker="*",
        label="DISK // Impactor // Silicate",
        alpha=0.6,
        cmap='viridis'
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.ax.set_title("Radius")
    ax.set_xlabel("Entropy")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Disk Entropy vs. Temperature (Radius Colored)")
    ax.grid()
    ax.legend()
    plt.savefig("temperature_vs_entropy_and_radius.png", format='png')
