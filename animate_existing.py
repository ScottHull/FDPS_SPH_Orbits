#!/usr/bin/env python3
from src.animate import animate

min_iteration = 0
max_iteration = 3000
increment = 1
fps = 30
path = "500_b073_new"
to_name = "3D_{}.mp4".format(path)

animate(
    start_time=min_iteration,
    end_time=max_iteration,
    interval=increment,
    path=path,
    fps=fps,
    filename=to_name,
)
