import glob
import numpy as np
import pathlib
import argparse

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
    
    # gets the firing rate
    def get_fr(self):
        print()

if __name__ == '__main__':
    args = parser.parse_args()

    mouse_id = args.mouseID
    probe = args.probe
    qc = qcChecker(mouse_id, probe)