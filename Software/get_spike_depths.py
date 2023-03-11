import numpy as np
from generate_metrics_paths import generate_metrics_path_ephys
import pathlib
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')
parser.add_argument('--probe', help='Desired probe')

def get_sort_dir(base_path: pathlib.Path, mouse_id: str, probe: str):
    paths = generate_metrics_path_ephys(base_path, mouse_id)
    probe_files = paths[int(probe[1])]

    probe_file = [f for f in probe_files if 'probe{}'.format(probe[0]) in f][0]
    sort_dir = pathlib.Path(probe_file).parent.absolute()

    return sort_dir

# gets the spike depth for the given mouse and probe
def generate_spike_depth(base_path: pathlib.Path, mouse_id: str, probe: str):
    sort_dir = get_sort_dir(base_path, mouse_id, probe)
    print(sort_dir)
    output_path = pathlib.Path(sort_dir, 'spike_depths.npy')
    sparse_features = np.load(os.path.join(sort_dir, 'pc_features.npy'), mmap_mode='r').squeeze().transpose((0, 2, 1))
    print(sparse_features.shape)
    cols = np.load(os.path.join(sort_dir, 'pc_feature_ind.npy'), mmap_mode='r')
    spike_templates = np.load(os.path.join(sort_dir, 'spike_templates.npy'), mmap_mode='r')
    spike_times = np.load(os.path.join(sort_dir, 'spike_times.npy'), mmap_mode='r')
    channel_positions = np.load(os.path.join(sort_dir, 'channel_positions.npy'), mmap_mode='r')

    nbatch = 50000
    c = 0
    spikes_depths = np.zeros_like(spike_times) * np.nan
    nspi = spikes_depths.shape[0]

    while True:
        ispi = np.arange(c, min(c + nbatch, nspi))
        # take only first component
        features = sparse_features[ispi, :, 0]
        features = np.maximum(features, 0) ** 2  # takes only positive values into account


        ichannels = cols[spike_templates[ispi]].astype(np.uint32)
        # features = np.square(self.sparse_features.data[ispi, :, 0])
        # ichannels = self.sparse_features.cols[self.spike_templates[ispi]].astype(np.int64)
        ypos = channel_positions[ichannels, 1]
        #ypos = ypos[:, 0, :]
        with np.errstate(divide='ignore'):
            print('Features', features.shape)
            print('Ypos', ypos.shape)
            spikes_depths[ispi, :] = (np.sum(np.transpose(ypos * features) /np.sum(features, axis=1), axis=0))
        c += nbatch
        if c >= nspi:
            break

    np.save(output_path, spikes_depths)

if __name__ == '__main__':
    args = parser.parse_args()

    mouse_id = args.mouseID
    probe = args.probe
    base_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    t_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')

    #output_path = pathlib.Path(t_path, mouse_id, '{}_kilo_data'.format(probe), 'spikes.depths.npy')
    generate_spike_depth(base_path, mouse_id, probe)