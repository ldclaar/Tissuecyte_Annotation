import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk
import argparse
import pathlib
import numpy as np
from pandas_profiling import ProfileReport

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)
parser.add_argument('--probe', help='Probe to be analyzed', required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    mouse_id = args.mouseID
    probe = args.probe
    """
    
    df_final_a1 = pd.read_csv('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe_A1_channels_626791_warped.csv'.format(mouse_id))
    reference = sitk.ReadImage('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/field_reference/average_template_25.nrrd')
    arr = sitk.GetArrayFromImage(reference).T

    print(df_final_a1)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(df_final_a1['ML'], df_final_a1['DV'], s=15 ,alpha=0.95)

    ax.imshow(arr.sum(axis=0), cmap='gray')

    plt.show()
    """
    """
    wave = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.ProbeA-AP/waveform_metrics.csv'))
    metrics = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.ProbeA-AP/metrics_test.csv'))

    wave_metrics = wave.merge(metrics, on='cluster_id')

    plt.plot(list(range(len(wave_metrics['d_prime']))), wave_metrics['d_prime'])
    plt.title('D prime')
    plt.show()
    
    feature_amp = wave_metrics['amplitude']
    feature_fire = metrics['firing_rate']
    feature_max_drift = metrics['max_drift']
    feature_spread = wave_metrics['spread']
    feature_repolar = wave_metrics['repolarization_slope']
    feature_velocity_above = wave_metrics['velocity_above']
    feature_velocity_below = wave_metrics['velocity_below']

    feature_amp = feature_amp.rolling(113, win_type='boxcar').mean().dropna()
    print('Amp size', len(feature_amp))
    feature_fire = feature_fire.rolling(113, win_type='boxcar').mean().dropna()

    print('Fire size', len(feature_fire))
    #feature_max_drift = feature_max_drift.rolling(230).mean().dropna()
    feature_spread = feature_spread.rolling(113, win_type='boxcar').mean().dropna()
    print('Spread size', len(feature_spread))
    feature_repolar = feature_repolar.rolling(113).mean().dropna()
    print('Repolar size', len(feature_repolar))

    feature_velocity_above = feature_velocity_above.rolling(14).mean().dropna()
    print('Velocity Above size', len(feature_velocity_above))
    feature_velocity_below = feature_velocity_below.rolling(16).mean().dropna()
    print('Velocity Below size', len(feature_velocity_below))

    fig, ax = plt.subplots(1, 6)

    ax[0].plot(list(range(len(feature_amp))), feature_amp)
    ax[0].set_title('Amplitude')

    ax[1].plot(list(range(len(feature_fire))), feature_fire)
    ax[1].set_title('Firing Rate')

    ax[2].plot(list(range(len(feature_spread))), feature_spread)
    ax[2].set_title('Spread')

    ax[3].plot(list(range(len(feature_repolar))), feature_repolar)
    ax[3].set_title('Repolarization_Slope')

    ax[4].plot(list(range(len(feature_velocity_above))), feature_velocity_above)
    ax[4].set_title('Velocity Above')

    ax[5].plot(list(range(len(feature_velocity_below))), feature_velocity_below)
    ax[5].set_title('Velocity Below')

    plt.show()
    """
    metrics = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.{}-AP/waveform_metrics.csv'.format(probe)))
    profile = ProfileReport(metrics)
    profile.to_file("output.html")