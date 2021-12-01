#!/usr/bin/env python3
import os
import shutil
from src.report import BuildReports

start_time = 0
end_time = 3000
interval = 1
number_processes = 200
from_path = "/home/theia/scotthull/200k/gi_new_eos"
to_path = "/home/theia/scotthull/200k/formatted_gi_new_eos"
eos_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"

# if os.path.exists(to_path):
#     shutil.rmtree(to_path)
os.mkdir(to_path)

r = BuildReports(
    to_dir=to_path,
    from_dir=from_path,
    start_time=start_time,
    end_time=end_time,
    number_processes=number_processes,
    eos_phase_path=eos_phase_path,
    interval=interval,
    accessory_path=None
)
r.make_reports(
    mp_pool_size=5,
    accessory_path=None
)
