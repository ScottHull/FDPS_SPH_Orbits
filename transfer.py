#!/usr/bin/env python3
import os
import paramiko
from scp import SCPClient
from multiprocessing.pool import ThreadPool as Pool

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
    scp.put(from_path + i, to_path + i)
    os.remove(from_path + i)

# define worker function before a Pool is instantiated
def worker(scp, full_from_path, full_to_path):
    try:
        transfer_file(scp, full_from_path, full_to_path)
    except:
        print('error with item')

while True:
    pool = Pool(pool_size)
    for i in os.listdir(from_path):
        if ".dat" in i:
            pool.apply_async(worker, (scp, from_path + i, to_path + i,))
    pool.close()
    pool.join()
