import os
import numpy as np
import pandas as pd


for i in os.listdir(os.getcwd()):
    if i.endswith(".csv"):
        df = pd.read_csv(i)
        positions = list(zip(df['x'], df['y'], df['z']))
        velocities = list(zip(df['vx'], df['vy'], df['vz']))
        mass = list(df['mass'])
        df['angular_momentum'] = np.array([k * np.linalg.norm(np.cross(i, j)) for i, j, k in zip(positions, velocities, mass)])
        df.to_csv(i)
