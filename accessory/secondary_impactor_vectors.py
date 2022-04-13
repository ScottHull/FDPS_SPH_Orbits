from math import sqrt, atan, pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt

square_scale = 5
earth_radius = 1
theia_radius = 0.3
theia_imp_y = -0.7
earth = plt.Circle((0, 0), earth_radius, color='r')

def get_x_impact_loc(y):
    """
    Equation of a circle: x^2 + y^2 + r^2
    If Theia impacts the Earth at y = 0.7, then:
    x = sqrt(r^2 - y^2)
    :return:
    """
    # point of impact with the Earth
    return sqrt(earth_radius ** 2 - y ** 2)

def get_impact_angle(y):
    x = get_x_impact_loc(y)
    imp_point = np.array([x, y])
    return atan(imp_point[1] / imp_point[0])

def get_theia_position(y_imp):
    """
        Equation of a circle: x^2 + y^2 + r^2
        If Theia impacts the Earth at y = 0.7, then:
        x = sqrt(r^2 - y^2)
        :return:
        """
    imp_angle = get_impact_angle(y_imp)
    hypotenuse = theia_radius + earth_radius
    x_offset = hypotenuse * cos(imp_angle)
    y_offset = hypotenuse * sin(imp_angle)
    return np.array([x_offset, y_offset])



def theia(y_imp):
    """
    Equation of a circle: x^2 + y^2 + r^2
    If Theia impacts the Earth at y = 0.7, then:
    x = sqrt(r^2 - y^2)
    :return:
    """
    theia_center = get_theia_position(y_imp=y_imp)
    return plt.Circle((theia_center[0], theia_center[1]), theia_radius, color='b')


fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.set_xlim(-square_scale, square_scale)
ax.set_ylim(-square_scale, square_scale)

ax.add_patch(earth)
ax.add_patch(theia(y_imp=theia_imp_y))

earth_com_to_theia_com = get_theia_position(theia_imp_y) - np.array([0, 0])
impact_angle = get_impact_angle(y=theia_imp_y)
theia_imp_x = get_x_impact_loc(theia_imp_y)
sub_x, sub_y = earth_radius * cos(impact_angle), earth_radius * sin(impact_angle)
surface_normal = earth_com_to_theia_com - np.array([sub_x, sub_y])
ax.plot([earth_com_to_theia_com[0], 0], [earth_com_to_theia_com[1], 0], c='w')


ax.grid()
ax.set_aspect('equal')

plt.show()
