import os
import shutil
from math import sqrt
import matplotlib.pyplot as plt

from src.animate import animate

from src.reverse_time import ReverseTime

# target_path = "/Users/scotthull/Desktop/input/tar.dat"
# impactor_path = "/Users/scotthull/Desktop/input/imp.dat"
target_path = "/home/shull4/drift_test/FDPS_SPH/input/tar.dat"
impactor_path = "/home/shull4/drift_test/FDPS_SPH/input/imp.dat"
output_path = "/scratch/shull4/reversed_outputs"

if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.mkdir(output_path)

rt = ReverseTime(
    target_file_path=target_path,
    impactor_file_path=impactor_path,
    impact_parameter=0.73,
    dt=-5,
    center_target=True
)

v_esc = sqrt(2.0 * rt.G * (rt.target_mass + rt.impactor_mass) / (rt.radius_target + rt.radius_impactor))

loop = 0
plotted = 0
while rt.distance < 3 * rt.radius_target:
    rt.reverse()
    if loop % 5 == 0:
        fig = rt.plot_current_position()
        plt.savefig(output_path + "/{}.png".format(plotted), format='png')
        plotted += 1
    loop += 1
# fig = rt.plot_current_position()
# plt.show()

print(
    "INITIAL SETUP:\n"
    "IMPACTOR-TARGET DISTANCE: {} Rad Tar ({} m)\n"
    "IMPACTOR-TARGET X DISTANCE: {} Rad Tar ({} m)\n"
    "IMPACTOR-TARGET Y DISTANCE: {} Rad Tar ({} m)\n"
    "IMPACTOR-TARGET Z DISTANCE: {} Rad Tar ({} m)\n"
    "INITIAL TARGET COORDS: {}\n"
    "INITIAL IMPACTOR COORDS: {}\n"
    "ESCAPE VELOCITY: {}\n"
    "INITIAL_TARGET_X_VELOCITY: {} v_esc\n"
    "INITIAL_TARGET_Y_VELOCITY: {} v_esc\n"
    "INITIAL_TARGET_Z_VELOCITY: {} v_esc\n"
    "INITIAL_IMPACTOR_X_VELOCITY: {} v_esc\n"
    "INITIAL_IMPACTOR_Y_VELOCITY: {} v_esc\n"
    "INITIAL_IMPACTOR_Z_VELOCITY: {} v_esc\n".format(
        rt.distance / rt.radius_target,
        rt.x_distance / rt.radius_target, rt.y_distance / rt.radius_target, rt.z_distance / rt.radius_target,
        rt.distance, rt.x_distance, rt.y_distance, rt.z_distance,
        rt.com_target, rt.com_impactor, v_esc,
        rt.v_target_x / v_esc, rt.v_target_y / v_esc, rt.v_target_z / v_esc,
        rt.v_impactor_x / v_esc, rt.v_impactor_y / v_esc, rt.v_impactor_z / v_esc,
    )
)

animate(
    start_time=0,
    end_time=plotted - 1,
    interval=1,
    path=output_path,
    fps=10,
    filename="reverse_time.mp4",
    reverse=True
)

rt.plot_velocity_history()
