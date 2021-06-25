from math import log10
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'size': 18}

matplotlib.rc('font', **font)

def big_O_scaler(desired_particles):
    # N = 2 * 10^6 simulates 1 second in 0.3 wall clock seconds
    # we want to run 20 simulated hours (7200 simulated seconds)
    scaler = 0.3 * 7200  # (0.3 t_wall / 1 t_sim) * (7200 t_sim) = 6 days wall clock time
    # what if we want to find out how long N = 10^7 takes?
    # i.e. 6 days / (2 * 10^6 * log10(2 * 10^6)) = x / (10^7 * log10(10^7))
    return desired_particles * log10(desired_particles) * (21600 / ((2 * 10 ** 6) * log10(2 * 10 ** 6)))

def seconds_to_days(s):
    return s / 86400


r = np.arange(10, 10 ** 8, 1000)
d = [big_O_scaler(desired_particles=i) for i in r]
d = [seconds_to_days(s=i) for i in d]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    r, d, linewidth=2.0, color='black'
)
ax.set_xlabel("N")
ax.set_ylabel("Days")
ax.set_title("Computational Time of FDPS SPH (O = N log10(N))")
ax.grid()

plt.show()
