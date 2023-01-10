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
    print(path)
    parent = path.parent.absolute()
    mean_waveforms = np.load(os.path.join(parent, 'mean_waveforms.npy'))
    amps = np.load(os.path.join(parent, 'amplitudes.npy'))
    print(amps.shape)
    print(mean_waveforms.shape)

    m_w_proj = mean_waveforms.sum(axis=0)
    print(m_w_proj.shape)

    fig, ax = plt.subplots(1, 1)

    ax.plot(list(range(len(m_w_proj))), m_w_proj)
    ax.set_title('Mean Waveforms')
    #ax[1].plot(list(range(len(amps))), amps)
    #ax[1].set_title('Raw Amplitudes')
    plt.show()
