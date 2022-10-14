# resamples the image series
# input - the location of the series
# input - the output directory of the resampled images

import SimpleITK as sitk
import os
import subprocess
from subprocess import Popen, PIPE
import argparse
import pathlib
from psycopg2 import connect, extras
import paramiko
import socket

parser = argparse.ArgumentParser()
parser.add_argument('--user', help='Allen Insitute User name', required=True)
parser.add_argument('--password', help='ALlen Institute Password', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
parser.add_argument('--downsampleFactor', help='Downsample images by this factor', required=False)

# affine aligns the image series using the input xml and transform xml
# first ssh into hpc using credentials provided
# allocate resources
# run affine aligned process
# saves the output to a location on the network
def affine_aligned(resampler, inputXMLFile, transformXMLFile, zeroFieldFile, outputDirectory) :
    cmd = [resampler,  '--downsampleFactor', '8', inputXMLFile, transformXMLFile, zeroFieldFile, outputDirectory]
    subprocess.call(cmd)

# ssh into hpc using username and password
# calls resample shell script to do the resampling with given mouse id
# add downsampling factor
def ssh_to_hpc(user, psswd, mouse_id):
    #proc = Popen(['ssh', '-T', 'arjun.sridhar@hpc-login'], stdin=PIPE, encoding='utf8')
    #proc.communicate(input='$otex5Xp')

    print('User', user)
    hostname='hpc-login'
    username=user
    password=psswd

    #cmd_allocation = 'srun -c 1 --mem=1mb -p celltypes --pty bash; pwd'
    # cd to directory and and call resampling job
    cmd_execute = 'cd /allen/scratch/aibstemp/arjun.sridhar; pwd; srun -N1 -c50 -t5:00:00 --mem=250gb -p celltypes resample_job.sh {}'.format(mouse_id) 

    ssh=paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname, username=username, password=password)
    print('Connection succesful, executing resampling')

    #stdin,stdout,stderr=ssh.exec_command(cmd_allocation)
    #print('Allocation output', stdout.read().decode())
    #print('Allocation Error', stderr.read().decode())

    stdin,stdout,stderr=ssh.exec_command(cmd_execute)
    print('Output', stdout.read().decode())
    print(stderr.read().decode())
    ssh.close()

def get_specimen_id_from_labtracks_id(labtracks_id):
    SPECIMEN_QRY = '''
        SELECT *
        FROM specimens sp
        WHERE sp.external_specimen_name=cast({} as character varying)
    '''

    mouse_info = query_lims(SPECIMEN_QRY.format(int(labtracks_id)))
    return mouse_info[0]['id']

# query lims and return result from given query
def query_lims(query_string):
    con = connect(
        dbname='lims2',
        user='limsreader',
        host='limsdb2',
        password='limsro',
        port=5432,
    )
    con.set_session(
        readonly=True, 
        autocommit=True,
    )
    cursor = con.cursor(
        cursor_factory=extras.RealDictCursor,
    )
    cursor.execute(query_string)
    result = cursor.fetchall()
    
    return result

# gets the tissuecyte info for the mouse id
# ssh into cluster using credentials given
# runs resampling and outputs to directory
def get_tc_info(mouse_id):
    TISSUECYTE_QRY = '''
        SELECT *
        FROM image_series im
        WHERE im.specimen_id = {}
    '''

    tc = dict(query_lims(TISSUECYTE_QRY.format(get_specimen_id_from_labtracks_id(int(mouse_id)))))
    storage_directory = tc['storage_directory']

    inputXMLFile = os.path.join(storage_directory, 'local_alignment', 'localalignment_input.xml')
    outputXMLFile = os.path.join(storage_directory, 'local_alignment', 'localalignment_transform_input.xml')

    return inputXMLFile, outputXMLFile
if __name__ == '__main__':
    """
    # get the input arguments
    args = parser.parse_args()
    input_series = pathlib.Path(args.inputSeries)
    output_dir = pathlib.Path(args.outputDir)

    # directories
    model_directory = '/allen/programs/celltypes/production/0378/informatics/model_september_2017/P56'
    resampler = '/shared/bioapps.x86_64/infoapps/lims2_modules/mouseconn/alignment/bin/idpProjectionResampleVolumeRGBModule'
    working_directory = '/allen/scratch/aibstemp/arjun.sridhar/ForCorbett/image_series_1183416716'

    # get deformation field and reference atlas
    field_file = os.path.join( working_directory, 'output_10', 'deformation_upsampled.nrrd')
    reference_file = os.path.join(model_directory, 'atlasVolume', 'average_template_10.nrrd')

    reference = sitk.ReadImage(reference_file)
    field = sitk.ReadImage(field_file)

    # make a zero warp field
    zeroField = sitk.Image(field.GetSize(), field.GetPixelID())
    zeroField.CopyInformation(field)

    output_file = os.path.join(working_directory, 'output_10', 'zeroField.mhd')
    sitk.WriteImage(zeroField, output_file, True)

    # create versions of resampled red, green, blue with no deformation
    inputXMLFile = os.path.join(input_series, 'grid', 'localalignment_input.xml' )
    outputXMLFile = os.path.join(input_series, 'grid','localalignment_transform_input.xml')
    output_directory = os.path.join(output_dir, 'resampled_images')

    affine_aligned(inputXMLFile, outputXMLFile, output_file, output_directory)
    """

    """
    # directories
    args = parser.parse_args()
    mouse_id = args.mouseID

    inputXMLFile, outputXMLFile = get_tc_info(mouse_id)
    model_directory = '/allen/programs/celltypes/production/0378/informatics/model_september_2017/P56'
    resampler = '/shared/bioapps.x86_64/infoapps/lims2_modules/mouseconn/alignment/bin/idpProjectionResampleVolumeRGBModule'
    working_directory = '/allen/scratch/aibstemp/arjun.sridhar/ForCorbett/image_series_1183416716'

    # get deformation field and reference atlas
    field_file = os.path.join( working_directory, 'output_10', 'deformation_upsampled.nrrd')
    reference_file = os.path.join(model_directory, 'atlasVolume', 'average_template_10.nrrd')

    reference = sitk.ReadImage(reference_file)
    field = sitk.ReadImage(field_file)

    # make a zero warp field
    zeroField = sitk.Image(field.GetSize(), field.GetPixelID())
    zeroField.CopyInformation(field)

    output_file = os.path.join(working_directory, 'output_10', 'zeroField.mhd')
    output_dir = '/allen/scratch/aibstemp/arjun.sridhar/resampled_temp'
    sitk.WriteImage(zeroField, output_file, True)
    affine_aligned(resampler, inputXMLFile, outputXMLFile, output_file, output_dir)
    """
    args = parser.parse_args()
    user = args.user
    password = args.password
    mouse_id = args.mouseID
    ssh_to_hpc(user, password, mouse_id)