import visvis as vis
import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk
import argparse
import pathlib
import numpy as np
import pandas as pd
#from pandas_profiling import ProfileReport
import os
import pickle
from generate_metrics_paths import generate_metrics_path_days
import pathlib
from scipy.spatial.transform import Rotation as R


def update_plot(linepts):
    vol = np.rot90(volume, k=2, axes=(0,2))
    vis.volshow3(vol)
    temp = volume.shape[2] 

    r = R.from_rotvec([0, np.pi, 0])
    linepts[:, 0] -= int(volume.shape[0] / 2)
    linepts[:, 1] -= int(volume.shape[1] / 2)
    linepts[:, 2] -= int(volume.shape[2] / 2)
    linepts = r.apply(linepts)

    linepts[:, 0] += int(volume.shape[0] / 2)
    linepts[:, 1] += int(volume.shape[1] / 2)
    linepts[:, 2] += int(volume.shape[2] / 2)

    vis.plot(temp - linepts[:, 2], linepts[:, 1], linepts[:, 0], lw=3, lc=(1, 0, 0))

if __name__ =='__main__':
    workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
    with open(os.path.join(workingDirectory, 'field_reference', 'acrnm_map.pkl'), 'rb') as f:
        acrnm_map = pickle.load(f)

    model_directory = pathlib.Path('{}/field_reference'.format(workingDirectory))
    reference_file = os.path.join(model_directory, 'average_template_25.nrrd')
    volume = sitk.GetArrayFromImage(sitk.ReadImage(reference_file)).T
    struct_acrnm = acrnm_map['VISpm2/3']
    anno = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(workingDirectory, 'field_reference', 'ccf_ano.mhd')))

    linepts = np.argwhere(anno == struct_acrnm)

    fig = vis.figure()

    # Create first axes
    a1 = vis.subplot(121)
    update_plot(linepts)
    # Enter main loop
    app = vis.use() # let visvis chose a backend for me
    app.Run()