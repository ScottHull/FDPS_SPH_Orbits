
def cross_section_xy(particles, min_z, max_z):
    """
    Want to hold z constant, between two values of Z.
    :param particles:
    :return:
    """
    return [
        p for p in particles if min_z <= p.position[2] <= max_z
    ]

def sort_particles_by_closest(particles):
    """
    Sorts particles in the z-direction, such that the particles with the greatest z value can be plotted first.
    """
    return [x for _, x in sorted(zip([p.position[2] for p in particles], particles))]

def grid_particles
