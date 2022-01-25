import os
import csv
import matplotlib.pyplot as plt

path = "/home/theia/scotthull/Paper1_SPH/gi/500_b073_new/formatted_500_b073_new"
output = open("times_vs_iteration.csv", 'w')
output.write("iteration,time\n")

def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours

iterations, times = [], []
for f in os.listdir(path):
    iteration = int(f.replace(".csv", ""))
    iterations.append(iteration)
    fname = path + "/{}".format(f)
    time = get_time(fname)
    times.append(time)
    output.write("{},{}\n".format(iteration, time))

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    iterations, times, linewidth=2.0
)
ax.grid()
ax.set_xlabel("Iteration")
ax.set_ylabel("Time (hrs)")
ax.set_title("Iteration vs. Time for All Simulations")
plt.savefig("iteration_vs_time.png", format='png')
