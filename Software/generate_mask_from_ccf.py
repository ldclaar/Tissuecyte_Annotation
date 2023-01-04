import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache, MouseConnectivityApi
from allensdk.api.queries.image_download_api import ImageDownloadApi
from PyQt5.QtCore import Qt, QAbstractTableModel
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, QTableView, \
    QVBoxLayout, QWidget, QPushButton, QGridLayout, QCheckBox, QFormLayout, QFileDialog, QComboBox, QMessageBox, QLineEdit
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtGui import QColor
import SimpleITK as sitk
import os
from get_tissuecyte_info import get_tc_info
import pathlib
import xmltodict
import argparse
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from PIL import Image
from sklearn.cluster import KMeans
from warp_image import warp_execute, cluster_annotations
import visvis as vis
import pickle

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session')
parser.add_argument('--probe', help='Probe to get labels for')

backend = 'pyqt5'
app = vis.use(backend)

DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]
class volumeMask(QWidget):
    def __init__(self, mouse_id, probe):
        super().__init__()
        self.mouseID = mouse_id
        # important directories
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.mcc = MouseConnectivityCache(resolution=25)
        self.anno, self.meta = self.mcc.get_annotation_volume()

        im = sitk.GetImageFromArray(self.anno)
        sitk.WriteImage(im, '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/field_reference/ccf_ano.mhd')
        self.rsp = self.mcc.get_reference_space()
        self.rsp.remove_unassigned(); # This removes ids that are not in this particular reference space
        self.name_map = self.rsp.structure_tree.get_name_map() # dictionary mapping ids to structure names
        self.acrnm_map = self.rsp.structure_tree.get_id_acronym_map() # dictionary mapping acronyms to ids
        self.colormap = self.rsp.structure_tree.get_colormap() # the colormap used for the allen 3D mouse atlas ontology

        with open(os.path.join(self.workingDirectory, 'field_reference', 'name_map.pkl'), 'wb') as f:
            pickle.dump(self.name_map, f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'acrnm_map.pkl'), 'wb') as f:
            pickle.dump(self.acrnm_map, f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'color_map.pkl'), 'wb') as f:
            pickle.dump(self.colormap, f)


        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_file)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')

        self.reference = sitk.ReadImage( self.reference_file )
        self.field = sitk.ReadImage( self.field_file )
        #self.dfmfld_transform = sitk.DisplacementFieldTransform(self.field)

        self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))
        self.volumeRed = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_red.mhd'))).T
        #self.volumeGreen = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_green.mhd'))).T
        #self.volumeBlue = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_blue.mhd'))).T
        Figure = app.GetFigureClass()
        self.fig = Figure(self)
        self.ax = vis.subplot(111)

        self.mainLayout = QHBoxLayout()
        self.mainLayout.addWidget(self.fig._widget)
        self.updateDisplay(probe.replace('_', ' '))
        #self.generateMask()

        self.setLayout(self.mainLayout)
        self.showMaximized()

    def generateMask(self):
        points = np.array(np.where(self.volArray)).T.astype(float)
        dfp = np.array(list(map(self.dfmfld_transform.TransformPoint, points*25)))
        sections = [int(np.round(v[0] / 10)) for v in dfp] # get slice?
        sections = sorted(sections)
        areas = set()

        model = KMeans(n_clusters=10)
        model.fit_predict(np.array(sections).reshape(-1, 1))
        
        for center in model.cluster_centers_:
            s = int(np.round(center))
            if s in self.name_map:
                key = list(self.acrnm_map.keys())[list(self.acrnm_map.values()).index(s)]
                print(self.name_map[s])

        """
        for s in slices:
            if s in self.name_map:
                key = list(self.acrnm_map.keys())[list(self.acrnm_map.values()).index(s)]
                areas.add(self.name_map[s])

        for area in areas:
            print(area)
        """

    def getColorVolume(self, rgb_levels=DEFAULT_COLOR_VALUES):
        level_adjusted_arrays = []
        for colori, int_level in zip(['red', 'green', 'blue'], rgb_levels):
            colarray = np.clip(self.int_arrays[colori], a_min=int_level[0], a_max=int_level[1]) - int_level[0]
            colarray = (colarray * 255. / (int_level[1] - int_level[0])).astype('uint8')
            level_adjusted_arrays.append(colarray)
        return np.stack(level_adjusted_arrays, axis=-1)

    def sliceCCF(self, nonzero_points):
        print()

    # displays the region along the probe track
    # probe: the probe to be displayed from the drop down
    def updateDisplay(self, probe):
        #print('Probe', probe)
        x = self.probeAnnotations[self.probeAnnotations.probe_name == probe].ML
        y = self.probeAnnotations[self.probeAnnotations.probe_name == probe].DV
        
        z = self.probeAnnotations[self.probeAnnotations.probe_name == probe].AP

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
        """
        intensity_values_red = np.zeros((linepts.shape[0],160))
        intensity_values_green = np.zeros((linepts.shape[0],160))
        intensity_values_blue = np.zeros((linepts.shape[0],160))
        """

        self.coords = {}
        self.linepts = linepts
        self.volumeWarp = np.zeros((528, 320, 456))
        #self.intensityValues = intensity_values_red
        pos_struct = {}

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2])

            for k in range(-80,80):
                try:
                    self.volumeWarp[int(np.round(linepts[j,0] / 2.5)),int(np.round(linepts[j,1]) / 2.5),int(np.round((linepts[j,2]+k) / 2.5))] = 1
                    #intensity_values_red[j,k+80] = (self.volumeRed[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    #intensity_values_green[j,k+80] = (self.volumeGreen[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    #intensity_values_blue[j,k+80] = (self.volumeBlue[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])

                except IndexError:
                    pass
        
        # display image
        #intensity_values *= 2.
        #self.trigScale = intensity_values.shape[0] / intensity_values.shape[1] 
        num_nonzero = len(np.nonzero(self.volumeWarp)[0])
        
        volume = sitk.GetImageFromArray(self.volumeWarp)
        volume.SetSpacing((25, 25, 25))
        dict_warp = {'AP': [], 'DV': [], 'ML': []}
        warped_image = warp_execute(volume, self.reference, self.field, sitk.sitkLinear)
        result = sitk.GetArrayFromImage(warped_image).T
        nonzero_points = np.argwhere(result > 0)

        structure_ids = set()
        struct_pos = {}
        for j in range(nonzero_points.shape[0]):
            try:
                structure_ids.add(self.anno[nonzero_points[j, 0], nonzero_points[j, 1], nonzero_points[j, 2]])
                struct_pos[self.anno[nonzero_points[j, 0], nonzero_points[j, 1], nonzero_points[j, 2]]] = (nonzero_points[j, 0], nonzero_points[j, 1], nonzero_points[j, 2])
            except IndexError:
                pass

        structs = []
        for struct in structure_ids:
            if struct in self.name_map:
                structs.append(self.name_map[struct])
        
        self.sliceCCF(nonzero_points)
        """
        vol3d = sitk.GetArrayFromImage(self.reference).T
        
        self.ax.bgcolor='k'
        self.ax.axis.axisColor = 'w'
        vis.cla()

        vis.volshow3(vol3d)
        legend = []
        values = list(self.acrnm_map.values())
        keys = list(self.acrnm_map.keys())

        for struct in structure_ids:
            if struct in self.colormap:
                col = tuple(self.colormap[struct])
                color = tuple(t / 255 for t in col)
                pos = struct_pos[struct]
                vis.plot(pos[2], pos[1], pos[0], mc=color, mw=5, ms='s', lw=0, mec=color, axes=self.ax)

                acrn = keys[values.index(struct)]
                legend.append(acrn)

        self.ax.legend = legend
        #self.ax.legend.bgcolor = 'grey'
        self.ax.axis.xLabel = 'ML'
        self.ax.axis.yLabel = 'DV'
        self.ax.axis.zLabel = 'AP'

        #cluster_annotations(num_nonzero, nonzero_points, dict_warp, None)
        
        self.int_arrays = {}
        self.int_arrays['red'] = intensity_values_red
        self.int_arrays['green'] = intensity_values_green
        self.int_arrays['blue'] = intensity_values_blue
        self.volArray = self.getColorVolume()
        self.volArray = self.volArray[250:, :, :]
        """

if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID
    probe = args.probe
    app.Create()
    mask = volumeMask(mouse_id, probe)
    app.Run()
