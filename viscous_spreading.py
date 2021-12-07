#!/usr/bin/env python3
import pandas as pd
import numpy as np
from math import sqrt
from statistics import mean
import matplotlib.pyplot as plt

new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
new_eos_silicate_table = ""
start_time = 0
end_time = 2500
increment = 100


mean_artificial_viscosity = 1.5


