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
        self.theia_client = self.createSSHClient(self.theia_server, self.theia_user, self.theia_pw)
        self.theia_scp = SCPClient(self.theia_client.get_transport())

    def createSSHClient(self, server, user, password):
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(server, username=user, password=password, port=22)
        return client

    def is_sftp_dir_exists(self, client, path):
        """
        Checks if a folder exists on the remote server.
        :param path:
        :return:
        """
        try:
            client.open_sftp().stat(path)
            return True
        except Exception:
            return False

    def create_sftp_dir(self, client, path):
        sftp = client.open_sftp()
        try:
            sftp.mkdir(path)
        except IOError as exc:
            if not self.is_sftp_dir_exists(client, path):
                raise exc

    def get_file(self, client, path, fname):
        f = path + "/" + fname
        return client.open_sftp().open(f)

    def get_df_from_theia(self, path, fname, skiprows=2):
        return pd.read_csv(self.get_file(self.theia_client, path, fname), skiprows=skiprows)

    def send_file_to_theia(self, from_path, to_path, filename):
        self.create_sftp_dir(self.theia_client, to_path)
        self.theia_scp.put("{}/{}".format(from_path, filename), "{}/{}".format(to_path, filename))
