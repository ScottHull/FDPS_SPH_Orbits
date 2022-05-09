import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14, })

new_hugoniot_path = "src/phase_data/forstSTS__hugoniot.txt"
old_hugoniot_path = "src/phase_data/duniteN__hugoniot.txt"

headers = [
    "Density", "Pressure", "Temp.", "Energy", "Sound Speed", "Entropy", "Shock Speed", "Part. Speed", "Phase"
]

new_hugoniot = pd.read_fwf(new_hugoniot_path, sep="\t", skiprows=1, header=None, names=headers)
old_hugoniot = pd.read_fwf(old_hugoniot_path, sep="\t", skiprows=1, header=None, names=headers)
print(new_hugoniot)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    1 / new_hugoniot["Density"], new_hugoniot["Pressure"] / 10 ** 9, linewidth=2.0, label="Stewart M-ANEOS Forsterite"
)
ax.plot(
    1 / old_hugoniot["Density"], old_hugoniot["Pressure"] / 10 ** 9,  linewidth=2.0, label="GADGET M-ANEOS Forsterite"
)

ax.set_xlabel(r" 1/ $\rho$ (m$^3$/kg)")
ax.set_ylabel(r"$P$ (GPa)")
ax.set_title("Hugoniot")
ax.grid(alpha=0.4)
ax.legend()

plt.show()
