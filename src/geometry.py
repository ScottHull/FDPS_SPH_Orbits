from math import atan, pi

from src.centering import center_of_mass


def get_impact_geometry(particles):
    target = [p for p in particles if p.tag <= 1]
    impactor = [p for p in particles if p.tag > 1]

    target_com_x, target_com_y, target_com_z = center_of_mass(
        x_coords=[p.position[0] for p in target],
        y_coords=[p.position[1] for p in target],
        z_coords=[p.position[2] for p in target],
        masses=[p.mass for p in target]
    )
    impactor_com_x, impactor_com_y, impactor_com_z = center_of_mass(
        x_coords=[p.position[0] for p in impactor],
        y_coords=[p.position[1] for p in impactor],
        z_coords=[p.position[2] for p in impactor],
        masses=[p.mass for p in impactor]
    )

    x_offset = impactor_com_x - target_com_x
    y_offset = impactor_com_y - target_com_y
    imp_angle = atan(y_offset / x_offset) * (180.0 / pi)

    min_x_tar = min([p.position[0] for p in target])
    max_x_tar = max([p.position[0] for p in target])

    min_y_tar = min([p.position[1] for p in target])
    max_y_tar = max([p.position[1] for p in target])

    min_x_imp = min([p.position[0] for p in impactor])
    max_x_imp = max([p.position[0] for p in impactor])

    min_y_imp = min([p.position[1] for p in impactor])
    max_y_imp = max([p.position[1] for p in impactor])

    return target, impactor, target_com_x, target_com_y, target_com_z, impactor_com_x, impactor_com_y, impactor_com_z, \
           imp_angle, min_x_tar, max_x_tar, min_y_tar, max_y_tar, min_x_imp, max_x_imp, min_y_imp, max_y_imp
