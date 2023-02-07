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
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)

if __name__ =='__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID
    basePath = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    
    waveMetricsPath = generate_metrics_path_days(basePath, mouse_id)
    days = sorted(list(waveMetricsPath.keys()))
    b2 = waveMetricsPath[days[1]][1]
    metrics = pd.read_csv(b2)
    plt.plot(metrics['velocity_above'], list(range(len(metrics['velocity_above']))))
    plt.show()

    parent = pathlib.Path(b2).parent.absolute()

    ksl = pd.read_table(os.path.join(parent, 'cluster_KSLabel.tsv'))
    print(ksl)

    channel_map = np.load(os.path.join(parent, 'channel_positions.npy'))
    print(channel_map)

    cluster_amp = pd.read_table(os.path.join(parent, 'cluster_Amplitude.tsv'))
    print(cluster_amp)

