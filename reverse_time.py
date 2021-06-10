import os
import shutil
import matplotlib.pyplot as plt

from src.animate import animate

from src.reverse_time import ReverseTime

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
    dt=-5
)

loop = 0
plotted = 0
while rt.distance < 3 * rt.radius_target:
    rt.reverse()
    if loop % 5 == 0:
        fig = rt.plot_current_position()
        plt.savefig(output_path + "/{}.png".format(plotted), format='png')
        plotted += 1
    loop += 1

print(
    "INITIAL SETUP:\n"
    "IMPACTOR-TARGET DISTANCE: {} Rad Tar\n"
    "IMPACTOR-TARGET X DISTANCE: {} Rad Tar\n"
    "IMPACTOR-TARGET Y DISTANCE: {} Rad Tar\n"
    "IMPACTOR-TARGET Z DISTANCE: {} Rad Tar\n"
    "INITIAL TARGET COORDS: {}\n"
    "INITIAL IMPACTOR COORDS: {}\n"
    "INITIAL_TARGET_X_VELOCITY: {}\n"
    "INITIAL_TARGET_Y_VELOCITY: {}\n"
    "INITIAL_IMPACTOR_X_VELOCITY: {}\n"
    "INITIAL_IMPACTOR_Y_VELOCITY: {}\n".format(
        rt.distance / rt.radius_target,
        rt.x_distance / rt.radius_target, rt.y_distance / rt.radius_target, rt.z_distance / rt.radius_target,
        rt.com_target, rt.com_impactor,
        rt.v_target_x, rt.v_target_y, rt.v_target_z,
        rt.v_impactor_x, rt.v_impactor_y, rt.v_impactor_z,
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

# rt.plot_velocity_history()
