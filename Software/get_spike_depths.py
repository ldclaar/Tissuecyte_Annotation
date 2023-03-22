import numpy as np
from generate_metrics_paths import generate_metrics_path_ephys, generate_metrics_path_days
import pathlib
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')
parser.add_argument('--probe', help='Desired probe')

def get_kilo_path_pilot(metrics_path: dict, probe_let_num: str):
    days = sorted(list(metrics_path.keys()))
    key = days[int(probe_let_num[1]) - 1]
    paths = metrics_path[key]
    path = ''

    for p in paths: # get the path for the metrics based on probe/day
        if 'probe' + probe_let_num[0] in p:
            path = p
            break

    return pathlib.Path(path).parent.absolute()

def get_templates(path:pathlib.Path):
    return np.load(pathlib.Path(path, 'templates.npy'), 'r+')

def get_channels(templates, cols):
    """ Gets peak channels for each waveform"""
    tmp = templates
    n_templates, n_samples, n_channels = tmp.shape
    if cols is None:
        template_peak_channels = np.argmax(tmp.max(axis=1) - tmp.min(axis=1), axis=1)
    else:
        # when the templates are sparse, the first channel is the highest amplitude channel
        template_peak_channels = cols[:, 0]

    return template_peak_channels

def get_templates_col(path:pathlib.Path):
    return np.load(pathlib.Path(path, 'templates_ind.npy'), 'r+')

def get_channel_positions(path:pathlib.Path):
    return np.load(pathlib.Path(path, 'channel_positions.npy'))

def get_spike_clusters(path:pathlib.Path):
    return np.load(pathlib.Path(path, 'spike_clusters.npy'), 'r+')

def generate_spike_depth_pilot(metrics_path, probe_let_num):
    path = get_kilo_path_pilot(metrics_path, probe_let_num)
    templates = get_templates(path)
    templates_cols = get_templates_col(path)
    templates_channels = get_channels(templates, templates_cols).astype(np.int64)
    channel_positions = get_channel_positions(path)
    spike_clusters = get_spike_clusters(path)
    
    clusters_depths = channel_positions[templates_channels, 1]
    spike_depths = clusters_depths[spike_clusters]
    np.save(pathlib.Path(path, 'spike_depths.npy'), spike_depths)

    return spike_depths

if __name__ == '__main__':
    args = parser.parse_args()

    mouse_id = args.mouseID
    probe = args.probe
    base_path = pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/datajoint/inbox/ks_paramset_idx_1')
    t_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')

    #output_path = pathlib.Path(t_path, mouse_id, '{}_kilo_data'.format(probe), 'spikes.depths.npy')
    metrics = generate_metrics_path_days(base_path, mouse_id)
    generate_spike_depth_pilot(metrics, probe)