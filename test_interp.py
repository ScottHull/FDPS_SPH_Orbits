import pandas as pd
from src.interpolation import GenericTrilinearInterpolation

path = "/Users/scotthull/Desktop/1000.csv"
eos = "src/phase_data/forst_STS.rho_u.txt"

eos_df = pd.read_fwf(eos, skiprows=2,
                     names=["density", "internal_energy", "temperature", "pressure", "soundspeed", "entropy"])
g = GenericTrilinearInterpolation(
    var1_array=list(eos_df['density']), var2_array=list(eos_df['internal_energy']),
    var3_array=list(eos_df['soundspeed'])
)

print(g.interpolate(var1=5000, var2=1e11))

