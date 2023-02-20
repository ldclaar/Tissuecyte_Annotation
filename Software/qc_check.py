import glob
import numpy as np
import pathlib
import argparse
from PyQt5 import QtGui
from matplotlib import cm

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session')
parser.add_argument('--probe', help='Desired probe')

# class for qc on kilosort output
class qcChecker():
    def __init__(self, mouse_id, probe):
        self.mouseID = mouse_id
        self.probe = probe
        self.kiloPath = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/{}_kilo_data'.format(mouse_id, probe))

        self.channelFiles = glob.glob(str(self.kiloPath) + '/channels*.npy')
        self.spikeFiles = glob.glob(str(self.kiloPath) + '/spikes*.npy')

        self.channelCoords = np.load(self.channelFiles[0]) # local coordinates file
        self.rawInd = np.load(self.channelFiles[1]) # raw ind file

        self.chnMin = np.min(self.channelCoords[:, 1])
        self.chnMax = np.max(self.channelCoords[:, 1])

        self.spikeAmps = np.load(self.spikeFiles[0]) # spikes amps file
        self.spikeClusters = np.load(self.spikeFiles[1]) # spike clusters file
        self.spikeDepths = np.load(self.spikeFiles[2]) # spike depths
        self.spikeTimes = np.load(self.spikeFiles[5]) # spike times
        self.spikeIdx = np.arange(self.spikeClusters.size)
        # Filter for nans in depths and also in amps
        self.kpIdx = np.where(~np.isnan(self.spikeDepths[self.spikeIdx]) &
                               ~np.isnan(self.spikeAmps[self.spikeIdx]))[0]
    
    # gets the amplitude
    def get_amp(self):
        A_BIN = 10
        amp_range = np.quantile(self.spikeAmps[self.spikeIdx][self.kpIdx], [0, 0.9])
        amp_bins = np.linspace(amp_range[0], amp_range[1], A_BIN)
        colour_bin = np.linspace(0.0, 1.0, A_BIN + 1)
        colours = ((cm.get_cmap('BuPu')(colour_bin)[np.newaxis, :, :3][0]) * 255).astype(np.int32)
        spikes_colours = np.empty(self.spikeAmps[self.spikeIdx][self.kpIdx].size,
                                      dtype=object)
        spikes_size = np.empty(self.spikeAmps[self.spikeIdx][self.kpIdx].size)
        for iA in range(amp_bins.size):
            if iA == (amp_bins.size - 1):
                idx = np.where((self.spikeAmps[self.spikeIdx][self.kpIdx] >
                                    amp_bins[iA]))[0]
                # Make saturated spikes a very dark purple
                spikes_colours[idx] = QtGui.QColor('#400080')
            else:
                idx = np.where((self.spikeAmps[self.spikeIdx][self.kpIdx] >
                                    amp_bins[iA]) &
                                   (self.spikeAmps[self.spikeIdx][self.kpIdx] <=
                                    amp_bins[iA + 1]))[0]
                spikes_colours[idx] = QtGui.QColor(*colours[iA])

            spikes_size[idx] = iA / (A_BIN / 4)

if __name__ == '__main__':
    args = parser.parse_args()

    mouse_id = args.mouseID
    probe = args.probe
    qc = qcChecker(mouse_id, probe)
    qc.get_amp()