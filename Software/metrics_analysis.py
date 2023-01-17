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
from sklearn.decomposition import PCA
from pandas_profiling import ProfileReport

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
#parser.add_argument('--probe', help='Probe to be analyzed', required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    mouse_id = args.mouseID

    basePath = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    waveMetricsPath = generate_metrics_path_days(basePath, mouse_id)
    days = list(sorted(waveMetricsPath.keys()))
    path = pathlib.Path(waveMetricsPath[days[1]][1])

    df = pd.read_csv(path)
    df = df.groupby('peak_channel').mean().reset_index()

    plt.plot(list(range(len(df['presence_ratio']))), df['presence_ratio'])
    plt.show()
    """
    profile = ProfileReport(df, explorative=True, minimal=True)
    profile.to_file('output.html')
    
    basePath = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    waveMetricsPath = generate_metrics_path_days(basePath, mouse_id)
    days = list(sorted(waveMetricsPath.keys()))
    path = pathlib.Path(waveMetricsPath[days[1]][1])
    print(path)
    parent = path.parent.absolute()
    mean_waveforms = np.load(os.path.join(parent, 'mean_waveforms.npy'))
    amps = np.load(os.path.join(parent, 'amplitudes.npy'))
    spike_clusters = np.load(os.path.join(parent, 'amplitudes.npy'))
    spike_times = np.load(os.path.join(parent, 'spike_times.npy'))


    print(amps.shape)
    print(mean_waveforms.shape)
    print(spike_clusters.shape)

    print(spike_clusters)


    fig, ax = plt.subplots(1, 2)

    ax[0].plot(list(range(len(amps))), amps)
    ax[0].set_title('Amps')
    ax[1].plot(list(range(len(spike_clusters))), spike_clusters)
    ax[1].set_title('Spike Clusters')
    plt.show()
    """
