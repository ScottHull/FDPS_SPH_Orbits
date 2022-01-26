#!/usr/bin/env python3
import os
import paramiko
from scp import SCPClient

class LunaToTheia:

    def __init__(self):
        self.theia_server = "epsl.earth.rochester.edu"
        self.theia_user = "scotthull"
        self.theia_pw = "Sunshine9249"
        # self.theia_client = self.createSSHClient(self.theia_server, self.theia_user, self.theia_pw).open_sftp()

    def createSSHClient(self, server, user, password):
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(server, username=user, password=password, port=22)
        return client

    def get_file(self, client, path, fname):
        f = path + "/" + fname
        return client.open(f)



