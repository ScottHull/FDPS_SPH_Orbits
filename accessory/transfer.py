#!/usr/bin/env python3
import os
import paramiko
from scp import SCPClient
import multiprocessing.dummy as mp

def createSSHClient(server, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, username=user, password=password, port=22)
    return client

base_folder_name = "gi_new_eos"
pool_size = 5  # your "parallelness"
from_path = "/scratch/shull4/{}/".format(base_folder_name)
to_path = "/home/theia/scotthull/{}/".format(base_folder_name)
server = "epsl.earth.rochester.edu"
user = "scotthull"
password = "PASSWORD"
ssh = createSSHClient(server, user, password)
scp = SCPClient(ssh.get_transport())

def transfer_file(scp, full_from_path, full_to_path):
    scp.put(full_from_path, full_to_path)
    os.remove(full_from_path)

# define worker function before a Pool is instantiated
def worker(args):
    inst_scp, full_from_path, full_to_path = args[0], args[1], args[2]
    try:
        transfer_file(inst_scp, full_from_path, full_to_path)
    except Exception as e:
        print(e)

while True:
    pool = mp.Pool(pool_size)
    pool.map(worker, [[scp, from_path + i, to_path + i] for i in os.listdir(from_path)])
    pool.close()
    pool.join()
