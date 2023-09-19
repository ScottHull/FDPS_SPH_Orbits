import numpy as np
import pandas as pd

from src.interpolation import NearestNeighbor1D


class GiantImpactReport:
    """
    This class is responsible for generating a report of the giant impact simulation.
    """

    def calculate_vmf(self, particles: pd.DataFrame, phase_curve: pd.DataFrame, nn: NearestNeighbor1D()):
        """
        Calculates the vapor mass fraction (vmf) of the particles.
        Assumes single phase particles.
        """
        # if the total entropy column is not defined, raise an error stating that it is required
        # ideally, total entropy should include both the particle entropy and the entropy gain associated with the
        # orbital circularization
        if 'total entropy' not in particles.columns:
            raise AttributeError("A 'total entropy' column is required to calculate the velocity mass fraction.\n"
                                 "This is the sum of the particle entropy and the entropy "
                                 "gain associated with orbital circularization.")
        # define the supercritical temperature
        supercritical_T = max(phase_curve['temperature'])
        # get the nearest temperature index for each particle on the phase curve corresponding to the
        # particle's temperature
        particles['nearest temperature index'] = nn.neighbor_index(given_val=particles['temperature'],
                                                                    array=list(phase_curve['temperature']))
        # get the corresponding phase curve entropy on the liquid and vapor curves
        particles['nearest liquid entropy'] = phase_curve['entropy_sol_liq'][particles['nearest temperature index']]
        particles['nearest vapor entropy'] = phase_curve['entropy_vap'][particles['nearest temperature index']]
        # define rules on whether the particle is pure liquid or pure vapor, or supercritical or mixed phase
        is_supercritical = particles['temperature'] >= supercritical_T
        is_pure_liquid = particles['total entropy'] <= particles['nearest liquid entropy']
        is_pure_vapor = particles['total entropy'] >= particles['nearest vapor entropy']
        # calculate the vmf for each particle
        particles['vmf'] = np.where(
            is_supercritical, 1.0, np.where(
                is_pure_liquid, 0.0, np.where(
                    is_pure_vapor, 1.0, (particles['total entropy'] - particles['nearest liquid entropy']) / (
                                particles['nearest vapor entropy'] - particles['nearest liquid entropy'])
                )
            )
        )
        return particles['vmf'].sum() / len(particles)
