from src.animate import animate

start = 0
end = 3000
interval = 5
fps = 10
from_path = "/home/theia/scotthull/FDPS_SPH_Orbits/new_and_old_animate"
fname = "animate.mp4"

animate(
    start_time=start,
    end_time=end,
    interval=interval,
    path=from_path,
    fps=fps,
    filename=fname,
)
