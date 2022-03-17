#!/usr/bin/env python3
import os
import shutil
from random import randint
import paramiko
from scp import SCPClient
import pandas as pd
import multiprocessing as mp

from src.combine import CombineFile


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
        sftp.close()

    def is_file_in_theia_path(self, path, fname):
        # can also do client.listdir()
        sftp = self.theia_client.open_sftp()
        if sftp.stat(path + "/" + fname):
            sftp.close()
            return True
        else:
            sftp.close()
            return False

    def listdir(self, client, path):
        sftp = client.open_sftp()
        l = sftp.listdir(path)
        sftp.close()
        return l

    def get_file(self, client, path, fname):
        f = path + "/" + fname
        return client.open_sftp().open(f)

    def get_df_from_theia(self, path, fname, skiprows=2):
        return pd.read_csv(self.get_file(self.theia_client, path, fname), skiprows=skiprows)

    def send_file_to_theia(self, from_path, to_path, filename):
        self.create_sftp_dir(self.theia_client, to_path)
        self.theia_scp.put("{}/{}".format(from_path, filename), "{}/{}".format(to_path, filename))

    def __get_single_output_filename(self, iteration, num_processes, current_process):
        file_format = "results.{}_{}_{}.dat"
        return file_format.format(str(iteration).zfill(5),
                                  str(num_processes).zfill(5),
                                  str(current_process).zfill(5))

    def get_file_from_theia(self, args):
        from_remote_path, to_local_path, fname, client = args
        ff = from_remote_path + "/{}".format(fname)
        ft = to_local_path + "/{}".format(fname)
        client.get(ff, ft)

    def get_and_combine_files_from_iteration(self, remote_path, num_processes, iteration,
                                             to_base_dir="/scratch/shull4"):
        ssh = self.createSSHClient(self.theia_server, self.theia_user, self.theia_pw)
        pool = mp.Pool(5)
        client = SCPClient(ssh.get_transport())
        to_path = to_base_dir + "/{}".format(randint(0, 100000))
        pool.map(self.get_file_from_theia, [[remote_path, to_path,
                                             self.__get_single_output_filename(iteration, num_processes, i),
                                             client] for i in range(0, num_processes)])
        pool.close()
        pool.join()

        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=num_processes, time=iteration, output_path=to_path, to_fname=to_fname)
        shutil.rmtree(to_path)
        return to_fname
