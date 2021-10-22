import paramiko
from scp import SCPClient

def createSSHClient(server, user, password, key_filename):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, user, password, key_filename=key_filename)
    return client

server = "epsl.earth.rochester.edu"
user = "scotthull"
password = "Sunshine9249"
key_filename = "id_rsa"
ssh = createSSHClient(server, user, password, key_filename=key_filename)
scp = SCPClient(ssh.get_transport())