#!/usr/bin/env python3
import os
import paramiko
from scp import SCPClient

def createSSHClient(server, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, username=user, password=password, port=22)
    return client

from_path = "/scratch/shull4/gi_new_eos/"
to_path = "/home/theia/scotthull/gi_new_eos/"
server = "epsl.earth.rochester.edu"
user = "scotthull"
password = "PASSWORD"
ssh = createSSHClient(server, user, password)
scp = SCPClient(ssh.get_transport())

while True:
    for i in os.listdir(from_path):
        if ".dat" in i:
            scp.put(from_path + i, to_path + i)
            os.remove(from_path + i)
