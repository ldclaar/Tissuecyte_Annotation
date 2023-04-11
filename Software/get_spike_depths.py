import numpy as np
import pathlib
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from qc_check import qcChecker

# gets the spike depth for a given kilo sort directory and saves numpy file in same directory
def generate_spike_depth(kilo_sort_dir:pathlib.Path) -> None:
    print(kilo_sort_dir)
    output_path = pathlib.Path(kilo_sort_dir, 'spike_depths.npy') # output path for spike depths
    sparse_features = np.load(os.path.join(kilo_sort_dir, 'pc_features.npy'), mmap_mode='r').squeeze().transpose((0, 2, 1)) # pc features
    print(sparse_features.shape)
    sparse_features_ind = np.load(os.path.join(kilo_sort_dir, 'pc_feature_ind.npy'), mmap_mode='r') # pc features ind
    spike_templates = np.load(os.path.join(kilo_sort_dir, 'spike_templates.npy'), mmap_mode='r')[:, 0]
    print(spike_templates.shape)
    spike_times = np.load(os.path.join(kilo_sort_dir, 'spike_times.npy'), mmap_mode='r')[:, 0]
    spike_times = spike_times / 30000.
    print(spike_times.shape)
    channel_positions = np.load(os.path.join(kilo_sort_dir, 'channel_positions.npy'), mmap_mode='r')

    nbatch = 50000
    c = 0
    spikes_depths = np.zeros_like(spike_times)
    nspi = spikes_depths.shape[0]

    while True:
        ispi = np.arange(c, min(c + nbatch, nspi))
        # take only first component
        features = sparse_features[ispi, :, 0]
        features = np.maximum(features, 0) ** 2  # takes only positive values into account


        ichannels = sparse_features_ind[spike_templates[ispi]].astype(np.uint32)
        # features = np.square(self.sparse_features.data[ispi, :, 0])
        # ichannels = self.sparse_features.cols[self.spike_templates[ispi]].astype(np.int64)
        ypos = channel_positions[ichannels, 1]
        #ypos = ypos[:, 0, :]
        with np.errstate(divide='ignore'):
            print('Features', features.shape)
            print('Ypos', ypos.shape)
            spikes_depths[ispi] = (np.sum(np.transpose(ypos * features) /np.sum(features, axis=1), axis=0))
        c += nbatch
        if c >= nspi:
            break

    np.save(output_path, spikes_depths)

if __name__ == '__main__':
    mouse_id = '644547'
    sessions = ['DRpilot_644866_20230207', 'DRpilot_644866_20230208', 'DRpilot_644866_20230209', 'DRpilot_644866_20230210']
    probes = ['probeA', 'probeB', 'probeC', 'probeD', 'probeE', 'probeF']

    for session in sessions:
        day = 1
        for probe in probes:
            path = '//allen/programs/mindscope/workgroups/templeton/TTOC/pilot recordings'.format(
                session, session, probe)
            kilo_sort_path = pathlib.Path(path)
            probe = probe[-1] + str(day)
            if kilo_sort_path.exists():
                generate_spike_depth(kilo_sort_path)
                qcChecker(kilo_sort_path, mouse_id, probe).get_correlation_data_img()

        day += 1