import numpy as np
import moviepy.editor as mpy


def animate(start_time, end_time, interval, path, filename="animation.mp4", fps=30, reverse=False):
    frames = [path + "/{}.png".format(time) for time in np.arange(start_time, end_time + interval, interval)]
    if reverse:
        frames = list(reversed(frames))
    animation = mpy.ImageSequenceClip(frames, fps=fps, load_images=True)
    animation.write_videofile(filename, fps=fps)


