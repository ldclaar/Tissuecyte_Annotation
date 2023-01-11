import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk
import argparse
import pathlib
import numpy as np
import pandas as pd
#from pandas_profiling import ProfileReport
import os
from generate_metrics_paths import generate_metrics_path_days
import pickle 

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
#parser.add_argument('--probe', help='Probe to be analyzed', required=True)

if __name__ == '__main__':
    workingDirectory = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte'
    anno = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(workingDirectory, 'field_reference', 'ccf_ano.mhd')))

    vol = np.zeros((528, 320, 456))
    with open(os.path.join(workingDirectory, 'field_reference', 'name_map.pkl'), 'rb') as f:
        name_map = pickle.load(f)
    with open(os.path.join(workingDirectory, 'field_reference', 'acrnm_map.pkl'), 'rb') as f:
        acrnm_map = pickle.load(f)
    with open(os.path.join(workingDirectory, 'field_reference', 'color_map.pkl'), 'rb') as f:
        colormap = pickle.load(f)
    df = pd.read_csv('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/608671/Probe_B2_channels_608671_warped.csv')
    channel_areas = []
    areas = set()
    keys = list(acrnm_map.keys())
    values = list(acrnm_map.values())
    for index, row in df.iterrows():
        struct = anno[row['AP'], row['DV'], row['ML']]

        if struct in values:
            ind = values.index(struct)
            key = keys[ind]
            channel_areas.append(key)
        else:
            channel_areas.append('N/A')
     
    df['channel_area'] = channel_areas
    print(df.loc[df['channel_area'] != 'N/A'])
