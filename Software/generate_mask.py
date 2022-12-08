import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache, MouseConnectivityApi
from allensdk.api.queries.image_download_api import ImageDownloadApi
import SimpleITK as sitk
import os
from get_tissuecyte_info import get_tc_info
import pathlib
import xmltodict
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session')

DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]
class volumeProperties:
    def __init__(self, mouse_id):
        self.mouseID = mouse_id
        self.dir = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte'
        self.workingDirectory = pathlib.Path('{}/{}'.format(self.dir, self.mouseID))
        print('Fetching Data')
        self.model_directory = pathlib.Path('{}/field_reference'.format(self.dir))
        self.probeAnnotations = pd.read_csv(os.path.join(self.workingDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))
        self.ccfImageDir = os.path.join(self.workingDirectory, '25_micron')

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))

        self.field = sitk.ReadImage( self.field_file, sitk.sitkVectorFloat64)

        self.mcc = MouseConnectivityCache(resolution=25)

        
        self.rsp = self.mcc.get_reference_space()
        self.rsp.remove_unassigned(); # This removes ids that are not in this particular reference space
        self.getAffineMatrix()
        self.name_map = self.rsp.structure_tree.get_name_map() # dictionary mapping ids to structure names
        self.acrnm_map = self.rsp.structure_tree.get_id_acronym_map() # dictionary mapping acronyms to ids
        self.colormap = self.rsp.structure_tree.get_colormap() # the colormap used for the allen 3D mouse atlas ontology
        self.regions_of_interest = self.get_structure_ids(list(self.acrnm_map.keys()))

        self.root_mask = self.rsp.make_structure_mask([997])
        self.root_points = np.array(np.where(self.root_mask)).T.astype(float)

        self.transformPoints()
        self.loadCCFVolume()
        self.updateDisplay('Probe B2')
        
    def loadCCFVolume(self):
        intensity_arrays = {}
        for imcolor in ['red']:
            resamp_image = sitk.ReadImage(os.path.join(self.ccfImageDir, 'resampled_{}.mhd'.format(imcolor)))
            intensity_arrays[imcolor] = sitk.GetArrayFromImage(resamp_image).T

        self.int_arrays = intensity_arrays

        self.volumeImage = self.getColorVolume()
        print('Data loaded')

    def transformPoints(self):
        # Transform root points
        dfmfld_transform = sitk.DisplacementFieldTransform(self.field)
        #aff_trans = sitk.AffineTransform(3)
        #aff_trans.SetMatrix(self.affine)

        dfp = np.array(list(map(dfmfld_transform.TransformPoint,self.root_points*25)))
        #self.volumeImage = dfp
        #tfp = np.array(list(map(aff_trans.TransformPoint,dfp)))

    # displays the region along the probe track
    # probe: the probe to be displayed from the drop down
    def updateDisplay(self, probe):
        x = self.probeAnnotations[self.probeAnnotations.probe_name == probe].ML / 2.5
        y = self.probeAnnotations[self.probeAnnotations.probe_name == probe].DV / 2.5
        
        z = self.probeAnnotations[self.probeAnnotations.probe_name == probe].AP / 2.5

        # get trajectory
        if len(z) > 0:
            data = np.vstack((z,y,x)).T
            datamean = data.mean(axis=0)

            D = data - datamean
            m1 = np.min(D[:,1]) * 2
            m2 = np.max(D[:,1]) * 2
            uu,dd,vv = np.linalg.svd(D)

            linepts = vv[0] * np.mgrid[-530:600:1][:,np.newaxis]
            linepts += datamean
            
            if linepts[-1,1] - linepts[0,1] < 0:
                linepts = np.flipud(linepts)
            
        # extract region of image, build row cells containing image
        intensity_values = np.zeros((linepts.shape[0],160))
        self.coords = {}
        self.linepts = linepts
        self.intensityValues = intensity_values

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2])
            #print(linepts[j, 0])
            for k in range(-80,80):
                try:
                    intensity_values[j,k+80] = (self.volumeImage[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                except IndexError:
                    pass
        
        # display image
        #intensity_values *= 2.
        #self.trigScale = intensity_values.shape[0] / intensity_values.shape[1] 
        #self.volArray = self.getColorVolume(intensity_values)
        self.volArray = intensity_values
        self.volArray = self.volArray[250:, :]

    def getColorVolume(self, rgb_levels=DEFAULT_COLOR_VALUES):
        level_adjusted_arrays = []
        for colori, int_level in zip(['red'], rgb_levels):
            colarray = np.clip(self.int_arrays[colori], a_min=int_level[0], a_max=int_level[1]) - int_level[0]
            colarray = (colarray * 255. / (int_level[1] - int_level[0])).astype('uint8')
            level_adjusted_arrays.append(colarray)
        return np.stack(level_adjusted_arrays, axis=-1)

    def get_structure_ids(self, acrnms):
        regions = []
        for acrnm in acrnms:
            regions.append(self.acrnm_map[acrnm])
        return regions

    def getAffineMatrix(self, direction='trv'):
        storage_directory = get_tc_info(self.mouseID)

        with open('/{}/local_alignment/localalignment_transform_input.xml'.format(storage_directory)) as f:
            transform_dict = xmltodict.parse(f.read())
            align3d = transform_dict['alignment']['alignment3ds']['alignment3d']
            #print(align3d)

            self.affine = np.zeros((4, 4), dtype=float)
            i = 0
            j = 0
            translation = False

            for key in align3d:
                if direction in key:
                    if not translation:
                        self.affine[i][j] = float(align3d[key]['#text'])
                        j += 1

                    if j == self.affine.shape[1] - 1 and not translation:
                        j = 0
                        i += 1

                    if i == self.affine.shape[0] - 1:
                        i = 0
                        translation = True
                        j = 3
                        continue

                    if translation:
                        self.affine[i][j] = float(align3d[key]['#text'])
                        i += 1

                    if i == self.affine.shape[0] - 1 and j == self.affine.shape[1] - 1:
                        break

            self.affine[-1][-1] = 1

if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID

    mask = volumeProperties(mouse_id)