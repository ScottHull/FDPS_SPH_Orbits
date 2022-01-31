import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("dark_background")

colors = print(plt.rcParams['axes.prop_cycle'])

N_eos = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH_Orbits/src/phase_data/duniteN__vapour_curve.txt"
STS_eos = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH_Orbits/src/phase_data/forstSTS__vapour_curve.txt"
f_cv_1_STS_eos = "/Users/scotthull/Desktop/forSTSM_vapour_curve.txt"  # f_cv = 1

N_eos_df = pd.read_fwf(N_eos, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

STS_eos_df = pd.read_fwf(STS_eos, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

f_cv_1_STS_eos_df = pd.read_fwf(f_cv_1_STS_eos, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    N_eos_df['entropy_sol_liq'],
    N_eos_df['temperature'],
    linewidth=2.0,
    # label="GADGET M-ANEOS Liquid",
    color='#8dd3c7'
)
ax.plot(
    N_eos_df['entropy_vap'],
    N_eos_df['temperature'],
    linewidth=2.0,
    # label="GADGET M-ANEOS Vapor",
    label="GADGET M-ANEOS",
    color='#8dd3c7'
)
ax.plot(
    STS_eos_df['entropy_sol_liq'],
    STS_eos_df['temperature'],
    linewidth=2.0,
    # label="Stewart M-ANEOS Liquid",
    color='#feffb3'
)
ax.plot(
    STS_eos_df['entropy_vap'],
    STS_eos_df['temperature'],
    linewidth=2.0,
    # label="Stewart M-ANEOS Vapor",
    label="Stewart M-ANEOS",
    color='#feffb3'
)
ax.plot(
    f_cv_1_STS_eos_df['entropy_sol_liq'],
    f_cv_1_STS_eos_df['temperature'],
    linewidth=2.0,
    # label="Stewart M-ANEOS Liquid ($f_{cv}$ = 1)",
    color='#bfbbd9'
)
ax.plot(
    f_cv_1_STS_eos_df['entropy_vap'],
    f_cv_1_STS_eos_df['temperature'],
    linewidth=2.0,
    # label="Stewart M-ANEOS Vapor ($f_{cv}$ = 1)",
    label="Stewart M-ANEOS ($f_{cv}$ = 1)",
    color='#bfbbd9'
)
ax.set_xlim(-50, 20000)
ax.set_xlabel("Entropy")
ax.set_ylabel("Temperature")
ax.set_title("M-ANEOS Boundary Comparision")
ax.grid(alpha=0.4)
ax.legend()
# plt.show()
plt.savefig("m-aneos_comparision.png", format='png', dpi=200)
