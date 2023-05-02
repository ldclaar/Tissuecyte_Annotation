# resamples the image series
# input - the location of the series
# input - the output directory of the resampled images

import SimpleITK as sitk
import os
import subprocess
from subprocess import Popen, PIPE
import argparse
import pathlib
import paramiko
import socket

parser = argparse.ArgumentParser()
parser.add_argument('--user', help='Allen Insitute User name', required=True)
parser.add_argument('--password', help='Allen Institute Password', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
parser.add_argument('--downsampleFactor', help='Downsample images by this factor', required=False)

# ssh into hpc using username and password
# calls resample shell script to do the resampling with given mouse id
# add downsampling factor
def ssh_to_hpc(user, psswd, mouse_id, down_sample):
    print('User', user)
    hostname='hpc-login'
    username=user
    password=psswd

    #cmd_allocation = 'srun -c 1 --mem=1mb -p celltypes --pty bash; pwd'
    # cd to directory and and call resampling job
    # cmd_execute = 'pwd; srun -N1 -c50 -t5:00:00 --mem=250gb -p celltypes /allen/programs/mindscope/workgroups/templeton-psychedelics/lesliec/resample_job.sh {} {}'.format(mouse_id, down_sample)
    cmd_execute = 'cd /allen/programs/mindscope/workgroups/templeton-psychedelics/lesliec; pwd; srun -N1 -c50 -t5:00:00 --mem=250gb -p celltypes resample_job.sh {} {}'.format(mouse_id, down_sample)

    ssh=paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname, username=username, password=password)
    print('Connection succesful, executing resampling')

    stdin,stdout,stderr=ssh.exec_command(cmd_execute)
    print('Output', stdout.read().decode())
    print(stderr.read().decode())
    ssh.close()

if __name__ == '__main__':
    args = parser.parse_args()
    user = args.user
    password = args.password
    mouse_id = args.mouseID

    if args.downsampleFactor is None:
        down_sample = '8'
    else:
        down_sample = args.downsampleFactor

    ssh_to_hpc(user, password, mouse_id, down_sample)
