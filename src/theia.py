#!/usr/bin/env python3
import os
import paramiko
from scp import SCPClient
import pandas as pd

class LunaToTheia:

    def __init__(self, server, u, p):
        self.theia_server = server
        self.theia_user = u
        self.theia_pw = p
        self.theia_client = self.createSSHClient(self.theia_server, self.theia_user, self.theia_pw).open_sftp()

    def createSSHClient(self, server, user, password):
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(server, username=user, password=password, port=22)
        return client

    def get_file(self, client, path, fname):
        f = path + "/" + fname
        return client.open(f)

    def get_df_from_theia(self, path, fname, skiprows=2):
        return pd.read_csv(self.get_file(self.theia_client, path, fname), skiprows=skiprows)



