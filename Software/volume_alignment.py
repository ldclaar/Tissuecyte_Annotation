from tkinter import N
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pathlib
import os
import SimpleITK as sitk
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, QProgressBar,\
    QVBoxLayout, QWidget, QPushButton, QGridLayout, QCheckBox, QFormLayout, QFileDialog, QComboBox, QMessageBox
from PyQt5 import QtCore
from pyqtgraph.Qt import QtGui
from PyQt5.QtGui import QColor
from PyQt5 import QtWidgets
import sys
from get_tissuecyte_info import get_tc_info
from warp_image import warp_channels
import xmltodict
from generate_metrics_paths import generate_metrics_path_days
import visvis as vis
from warp_image import warp_execute, cluster_annotations
import pickle
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)

SCALING_FACTOR = 1.5
DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]

# graph represents a metric plot that can be clicked on
# pyqtgraph documentation
class Graph(pg.GraphItem):
    def __init__(self):
        self.dragPoint = None
        self.dragOffset = None
        self.textItems = []
        self.pointPos = {}
        self.clicked = False
        pg.GraphItem.__init__(self)
        self.scatter.sigClicked.connect(self.onclick)

        
    def setData(self, **kwds):
        self.text = kwds.pop('text', [])
        self.data = kwds
        if 'pos' in self.data:
            npts = self.data['pos'].shape[0]
            self.data['data'] = np.empty(npts, dtype=[('index', int)])
            self.data['data']['index'] = np.arange(npts)
        self.setTexts(self.text)
        self.updateGraph()
        
    def setTexts(self, text):
        for i in self.textItems:
            i.scene().removeItem(i)
        self.textItems = []
        for t in text:
            item = pg.TextItem(t)
            self.textItems.append(item)
            item.setParentItem(self)
        
    def updateGraph(self):
        pg.GraphItem.setData(self, **self.data)
        for i,item in enumerate(self.textItems):
            item.setPos(*self.data['pos'][i])
       
    """
    def mouseDragEvent(self, ev):
        if ev.button() != QtCore.Qt.LeftButton:
            ev.ignore()
            return
        
        if ev.isStart():
            # We are already one step into the drag.
            # Find the point(s) at the mouse cursor when the button was first 
            # pressed:
            pos = ev.buttonDownPos()
            print('Pos', pos)
            pts = self.scatter.pointsAt(pos)
            if len(pts) == 0:
                ev.ignore()
                return
            self.dragPoint = pts[0]
            ind = pts[0].data()[0]
            self.dragOffset = self.data['pos'][ind] - pos
        elif ev.isFinish():
            self.dragPoint = None
            return
        else:
            if self.dragPoint is None:
                ev.ignore()
                return
        
        ind = self.dragPoint.data()[0]
        print(ind)
        self.data['pos'][ind] = ev.pos() + self.dragOffset
        self.updateGraph()
        ev.accept()
        self.pointPos[ind] = ev.pos()
    
    def clicked(self, pts):
        print("clicked: %s" % pts)
    """
    def onclick(self, plot, points):
        if hasattr(self, 'lastPoint'):
            self.lastPoint.setBrush(QtGui.QBrush(QColor('grey')))
            self.lastPoint.setPen(QtGui.QPen(QColor('grey')))

        self.clicked = True
        points[0].setBrush(QtGui.QBrush(QColor('green')))
        points[0].setPen(QtGui.QPen(QColor('green')))
        self.lastPoint = points[0]
        self.scatterPoint = points[0].pos()

# class for metric plot
# each metric plot is represented by a plot display item which uses the graph class above
class PlotDisplayItem():
    # create default image
    def __init__(self, measurement, waveform, metrics, probe_annotations, mouse_id, metrics_list):
        self.width = int(4000) #default = coronal view
        self.height = int(4000)
        self.remove = True
        self.show = True

        self.probeAnnotations = probe_annotations
        self.mouseID = mouse_id
        self.waveform = waveform
        self.metrics = metrics
        self.metricsList = metrics_list
        self.otherPlots = []
        self.oldChannels = [] # stack for undoing lines

        self.waveform_metrics = self.waveform.merge(self.metrics, on='cluster_id')
        self.measurement = measurement
        #self.generateMetrics(measurement)

        self.channelsPlot = Graph()
        #self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
        self.textItem = pg.TextItem(measurement.upper(), anchor=(1, 1))

        #self.channels = self.channelsOriginal
    
    # Updates the metrics with the given path to the metrics csv
    # path: string, the path to the csv
    def updateMetrics(self, path):
        self.waveform_metrics = pd.read_csv(path)
        print('Metrics Path', path)
        self.generateMetrics(self.measurement)

    # generates the metrics based on the measurement passed in
    # measurement: string, the metric
    def generateMetrics(self, measurement):
        if self.measurement == 'unit_density':
            self.processMetrics()
            self.frequencyCounts = self.waveform_metrics['peak_channel'].value_counts().sort_index().reset_index().to_numpy()
            #x = np.linspace(-10, 100, num=384)
            self.channelsRecorded = self.frequencyCounts[:, 0].tolist()
            self.channelsOriginal = []

            for i in range(384):
                if i in self.channelsRecorded:
                    ind = self.channelsRecorded.index(i)
                    self.channelsOriginal.append([i, self.frequencyCounts[ind, 1]])
                else:
                    self.channelsOriginal.append([i, 0])
           
            self.channelsOriginal = [[self.channelsOriginal[i][1], 384 - i - 1 + 256] for i in range(len(self.channelsOriginal))]
            
            x_val = [p[0] * 10 for p in self.channelsOriginal]
            conv = np.ones(10)

            smoothed = np.convolve(x_val, conv, mode='same')
            smoothed = smoothed / np.sum(conv)
            self.channelsOriginal = [[smoothed[i] - 100, self.channelsOriginal[i][1]] for i in range(384)]
        elif self.measurement == 'spread':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=0, shift_value=300)
        elif self.measurement == 'firing_rate':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=1/2, shift_value=250)
        elif self.measurement == 'd_prime':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=0, shift_value=250)
        elif self.measurement == 'cumulative_drift':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=0, shift_value=250)
        elif self.measurement == 'velocity_above':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=1/15, shift_value=250)
        elif self.measurement == 'velocity_below':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=1/15, shift_value=250)
        elif self.measurement == 'amplitude':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=5, shift_value=300)
        elif self.measurement == 'isi_viol':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=0, shift_value=250)
        elif self.measurement == 'nn_hit_rate':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=1/5, shift_value=250)
        elif self.measurement == 'presence_ratio':
            self.processMetrics()
            self.generateMetricChannels(measurement, scale_value=1/15, shift_value=250)

    # helper function to generate the metric channels
    # metric: string
    # scale_value: float, value to scale the plot by this amount when displayed
    # shit_value: int, value to shift the plot by this amount when displayed
    def generateMetricChannels(self, metric, scale_value, shift_value):
        #print('Metric', metric)
        peak_values = self.averageMetricsChannels['peak_channel'].values.tolist()
        values = self.averageMetricsChannels[metric]

        if 'velocity' in metric:
            values = values / scale_value
            #values = values / values.abs().max()
            values = values.to_numpy().tolist()
            conv = np.ones(10)
        else:
            conv = np.ones(10)
        #peak_values = [int(p) for p in peak_values]
        self.channelsOriginal = []

        for i in range(384):
            if i in peak_values:
                index = peak_values.index(i)
                if not pd.isna(values[index]):
                    self.channelsOriginal.append([values[index], 384 - i - 1 + 256])
                else:
                    self.channelsOriginal.append([0, 384 - i - 1 + 256])
                #conv[index] = 1
            else:
                self.channelsOriginal.append([0, 384 - i - 1 + 256])
        
        x_val = [p[0] for p in self.channelsOriginal]
        smoothed = np.convolve(x_val, conv, mode='same')
        smoothed = smoothed / np.sum(conv)
        #print(smoothed.shape)
        for i in range(384):
            if scale_value != 0 and 'velocity' not in metric:
                self.channelsOriginal[i] = [(smoothed[i] / scale_value) - shift_value, self.channelsOriginal[i][1]]
            else:
                self.channelsOriginal[i] = [(smoothed[i]) - shift_value, self.channelsOriginal[i][1]]
        
    def processMetrics(self):
        self.averageMetricsChannels = self.waveform_metrics.groupby('peak_channel').mean().reset_index()
        #self.averageMetricsChannels = self.averageMetricsChannels.rolling(rolling_value, win_type='boxcar').mean()

    # when a new probe is displayed
    def resetPlot(self, remove_probe=False):
        self.channels = self.channelsOriginal
        self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

    # helper function for shfiting points when new alignment is added
    # returns closest value to K
    # lst: list to look for closest
    # K: closest number to look for
    def closest(self, lst, K):
        if len(lst) > 0:
            return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

    # helper function to update channel positions for linear interpolation between 2 alignments
    # lp: list of points 
    # points_between: points that need to be updated with points from lp, same length as lp
    def replaceValues(self, lp, points_between):
        #print(points_between)
        for i in range(len(lp)):
            ind = self.channels.index(points_between[len(lp) - i - 1])
            #print(self.channels[ind])
            self.channels[ind][1] = lp[i]
            #print(self.channels[ind])

    # linearly space points between 2 alignments
    # pointsAdded: the list of alignments that have been made so far
    def linearSpacePoints(self, pointsAdded):
        if len(pointsAdded) > 1:
            sorted_linepts = sorted(pointsAdded)

            for i in range(len(sorted_linepts) - 1): # each time new line is added
                anchor_top = sorted_linepts[i]
                anchor_bottom = sorted_linepts[i + 1]

                #print(anchor_top, anchor_bottom)

                points_between = [p for p in self.channels if p[1] > anchor_top and p[1] < anchor_bottom] # current points between
                #print(len(points_between))
                lp = np.linspace(anchor_top + 1, anchor_bottom - 1, num=len(points_between)).tolist() # new linear interpolation
                #print(len(lp))
                #print(lp)
                #lp = [[p[0], int(p[1])] for p in lp]
                self.replaceValues(lp, points_between) # replace values with new interpolation
                #print(len(self.channels))

            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

    # for alignments above current channel
    # shift up, i.e. subtract
    # srt: list with sorted alignments already done
    # channel: y coord of point clicked on
    # diff: difference between channel and alignment 
    # flag: boolean, returns False if interpolation is done
    def linspaceAbove(self, srt, channel, diff, flag):
        for point in srt:
            if self.channelsPlot.scatterPoint[1] > point:
                points_between = [p for p in self.channels if p[1] < self.channelsPlot.scatterPoint[1] and p[1] > point]
                self.channels[channel][1] -= diff
                lp = np.linspace(point + 1, self.channels[channel][1] - 1, len(points_between))
                self.replaceValues(lp, points_between)
                flag = False
                break

        return flag

    # for alignments below current channel
    # shift down, i.e. add
    # srt: list with sorted alignments already done
    # channel: y coord of point clicked on
    # diff: difference between channel and alignment 
    # flag: boolean, returns False if interpolation is done
    def linspaceBelow(self, srt, channel, diff, flag):
        for point in srt:
            if self.channelsPlot.scatterPoint[1] < point:
                points_between = [p for p in self.channels if p[1] > self.channelsPlot.scatterPoint[1] and p[1] < point]
                self.channels[channel][1] += diff
                lp = np.linspace(self.channels[channel][1] + 1, point - 1, len(points_between))
                self.replaceValues(lp, points_between)
                flag = False
                break

        return flag

    # scales the plots using the given scale factors below
    def updateDisplay(self, probe, linepts, intensity_values, keep_y=False, old_y=None, points_added=None):
        if not keep_y:
            ap_scale = 1
            dv_scale = 1/0.94
            lr_scale = 1
 
            # get vectors in xyz
            self.dz = (linepts[-1, 0] - linepts[0, 0])
            self.dy = (linepts[-1, 1] - linepts[0, 1])
            self.dx = (linepts[-1, 2] - linepts[0, 2])
            self.vector = np.array([self.dz, self.dy, self.dx])
            
            # apply scale
            self.dzScale = self.dz * ap_scale
            self.dyScale = self.dy * dv_scale
            self.dxScale = self.dx * lr_scale
            self.vectorScale = np.array([self.dzScale, self.dyScale, self.dxScale])

            # get scale factor
            self.scale = np.linalg.norm(self.vectorScale) / np.linalg.norm(self.vector)
            #print(self.scale)
            #channels_scale = np.linspace(0, self.scale, 384)

            newPoints = []
            for i in range(len(self.channelsOriginal)):
                newPoints.append([self.channelsOriginal[i][0], self.channelsOriginal[i][1] * self.scale]) # apply scale factor
        
            #print('Initial Probe length in pixels', newPoints[0][1] - newPoints[-1][1])
            self.ogPoints = newPoints
            self.channels = newPoints
            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        else: # some alignment has been done already, so use existing coordinates
            newPoints = [] 
            for i in range(len(self.channelsOriginal)):
                newPoints.append([self.channelsOriginal[i][0], old_y[i][1]])

            self.ogPoints = newPoints
            self.channels = newPoints
            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
            self.linearSpacePoints(points_added)

# class to do the alignment between channels and regions
class VolumeAlignment(QWidget):
    def __init__(self, mouse_id):
        super().__init__()
        self.title = 'Volume Alignment'
        self.mouseID = mouse_id
        self.probe = 'ProbeA'
        self.waveform = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.{}-AP/waveform_metrics.csv'.format(self.probe)))
        self.metrics = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.{}-AP/metrics_test.csv'.format(self.probe)))
        
        self.basePath = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
        self.waveMetricsPath = generate_metrics_path_days(self.basePath, self.mouseID)
        self.days = sorted(list(self.waveMetricsPath.keys()))

        self.initUI()
        self.displayRegion()
    
    def initUI(self):
        self.setWindowTitle(self.title)
        self.width = int(4000) #default = coronal view
        self.height = int(4000)
        self.remove = True
        self.showProbe = True
        self.showMask = True
        self.prevProbe = ''

        im8 = np.ones((self.height,self.width),dtype='uint8')*255
        self.image = pg.image(im8)
        self.image.ui.histogram.hide()
        self.image.ui.roiBtn.hide()
        self.image.ui.menuBtn.hide()
        self.image.setObjectName('image')

        self.imageMask = pg.image(im8)
        self.imageMask.ui.histogram.hide()
        self.imageMask.ui.roiBtn.hide()
        self.imageMask.ui.menuBtn.hide()
        self.imageMask.setObjectName('image')
        #self.image.setImage(im8.transpose())

        self.pointsAdded = [] # y coords of alignments done
        self.channelCoords = {}
        self.lineItems = [] # pyqtgraph item for alignments
        self.oldChannels = [] # stack for undoing lines
        self.distPoints = []
        self.distItems = []
        self.s = []
        self.pos = []
        self.myFont = ImageFont.load_default()
        self.alignments = {}

        # important directories
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_file)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')
        
        
        with open(os.path.join(self.workingDirectory, 'field_reference', 'name_map.pkl'), 'rb') as f:
            self.name_map = pickle.load(f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'acrnm_map.pkl'), 'rb') as f:
            self.acrnm_map = pickle.load(f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'color_map.pkl'), 'rb') as f:
            self.colormap = pickle.load(f)
        self.anno = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.workingDirectory, 'field_reference', 'ccf_ano.mhd')))
        
        # read reference image, deformation field, probe annotations
        self.reference = sitk.ReadImage( self.reference_file )
        self.field = sitk.ReadImage( self.field_file)
        self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))

        # metrics list and plots to show
        self.metricsList = ['Spread', 'Amplitude', 'Isi_Viol', 'NN_Hit_Rate', 'Firing_Rate', 'Presence_Ratio', 'Velocity_Above', 'Velocity_Below']

        self.unitPlot = PlotDisplayItem('unit_density', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.spreadPlot = PlotDisplayItem('spread', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.isiPlot = PlotDisplayItem('isi_viol', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.firePlot = PlotDisplayItem('firing_rate', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.velocityAbovePlot = PlotDisplayItem('velocity_above', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.velocityBelowPlot = PlotDisplayItem('velocity_below', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.ampPlot = PlotDisplayItem('amplitude', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.hitRatePlot = PlotDisplayItem('nn_hit_rate', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.presPlot = PlotDisplayItem('presence_ratio', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        #self.repolarPlot = PlotDisplayItem('repolarization_slope', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID)
        #self.dPrimePlot = PlotDisplayItem('d_prime', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID)
        #self.cumDriftPlot = PlotDisplayItem('cumulative_drift', self.waveform, self.metrics, self.probeAnnotations, self.mouseID, self.metricsList)
        self.plots = {'unit_density': self.unitPlot, 'spread': self.spreadPlot, 'velocity_above': self.velocityAbovePlot, 'firing_rate': self.firePlot,
                      'velocity_below': self.velocityBelowPlot, 'amplitude': self.ampPlot, 'isi_viol': self.isiPlot, 'nn_hit_rate': self.hitRatePlot,
                      'presence_ratio': self.presPlot}
        
        for plot in self.plots:
            others = [k for k in self.plots if k != plot]
            for other in others:
                self.plots[plot].otherPlots.append(self.plots[other]) # add other plots to list for updating 

        #self.getColorVolume()

        # main layout with region of interest
        self.mainLayout = QVBoxLayout()
        self.imageLayout = QHBoxLayout()

        
        self.imageLayout.addWidget(self.image)
        self.imageLayout.addWidget(self.imageMask)

        # add ui features: probe/metric drop downs, red/green toggle, toggle probe, toggle mask, warp to ccf
        self.probes = self.probeAnnotations['probe_name'].unique()
        self.probeLetters = [s[s.index(' ')+1:][0] for s in self.probes]
        self.probeLetters = sorted(list(set(self.probeLetters)))
        print(self.probeLetters)
        self.probeDropDown = QComboBox()
        for probe in sorted(self.probes):
            self.probeDropDown.addItem(probe)

        self.probeDropDown.currentTextChanged.connect(self.displayRegion)
        self.redCheck = QCheckBox('Toggle Red')
        self.redCheck.setFocusPolicy(QtCore.Qt.NoFocus)
        self.redOld = None
        self.isRedChecked = False
        self.redCheck.clicked.connect(self.toggleRed)

        self.greenCheck = QCheckBox('Toggle Green')
        self.greenCheck.setFocusPolicy(QtCore.Qt.NoFocus)
        self.greenOld = None
        self.isGreenChecked = False
        self.greenCheck.clicked.connect(self.toggleGreen)

        self.metrics = QComboBox()
        self.metrics.currentTextChanged.connect(self.metricChanged)
        for metric in self.metricsList:
            self.metrics.addItem(metric)

        self.toggleProbeButton = QPushButton('Toggle Probe')
        self.toggleProbeButton.clicked.connect(self.toggleProbe)

        self.resetPlotButton = QPushButton('Reset Metric Plot')
        self.resetPlotButton.clicked.connect(self.resetPlot)

        self.viewButton = QPushButton('View Probe with Selected Metric')
        self.viewButton.clicked.connect(self.displayRegion)

        self.calcDistanceButton = QPushButton('Calculate Probe Distance')
        self.calcDistanceButton.clicked.connect(self.calcDistance)

        self.viewCCFButton = QPushButton('Toggle Mask')
        self.viewCCFButton.clicked.connect(self.toggleCCFRegions)
        self.showCCF = False

        self.warpButton = QPushButton('Warp to CCF')
        self.warpButton.clicked.connect(self.warpChannels)

        self.probeViewLayout = QHBoxLayout()
        self.probeViewLayout.addWidget(self.probeDropDown)
        self.probeViewLayout.addWidget(self.metrics)
        #self.probeViewLayout.addWidget(self.toggleProbeButton)
        self.probeViewLayout.addWidget(self.viewButton)
        self.probeViewLayout.addWidget(self.resetPlotButton)
        self.probeViewLayout.addWidget(self.toggleProbeButton)
        #self.probeViewLayout.addWidget(self.calcDistanceButton)
        self.probeViewLayout.addWidget(self.viewCCFButton)
        self.probeViewLayout.addWidget(self.redCheck)
        self.probeViewLayout.addWidget(self.greenCheck)

        self.probeViewLayout.addWidget(self.warpButton)
        self.probeViewLayout.setAlignment(QtCore.Qt.AlignTop)
        self.mainLayout.addLayout(self.imageLayout)
        self.mainLayout.addLayout(self.probeViewLayout)

        self.setLayout(self.mainLayout)
        self.showMaximized()

    # function to bring mask back 
    def showMaskHelper(self):
        im = Image.fromarray(self.volArray)
        overlay = Image.blend(im, self.mask, 0.5)
        #overlay.show()

        draw = ImageDraw.Draw(overlay)

        for label in self.labels_pos:
            positions = self.labels_pos[label]
            centroid = np.array(positions).mean(axis=0)
            draw.text((int(np.round(centroid[1])), int(np.round(centroid[0]))), label, font=self.myFont, fill=(255, 255, 255))

        self.blended = np.array(overlay)
        rot = np.rot90(self.blended)
        flip = np.flipud(rot)
        self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)

    # hides/shows the red channel of the volume
    def toggleRed(self):
        if not self.isRedChecked: # hide red channel
            self.redOld = self.volArray[:, :, 0].copy()
            self.volArray[:, :, 0] = 0
            self.isRedChecked = True

            if not self.showMask: # mask is currently hidden
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else: # bring back mask on to volume slice image
                self.showMaskHelper()
        else: # show red channel
            self.volArray[:, :, 0] = self.redOld.copy()
            self.isRedChecked = False

            if not self.showMask:
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else:
                self.showMaskHelper()

    # hides/shows the green channel of the volume
    def toggleGreen(self):
        if not self.isGreenChecked: # hide green
            self.greenOld = self.volArray[:, :, 1].copy()
            self.volArray[:, :, 1] = 0
            self.isGreenChecked = True

            if not self.showMask: # mask is currently hidden
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else: # bring back mask on to volume slice image
                self.showMaskHelper()
        else:
            self.volArray[:, :, 1] = self.greenOld.copy()
            self.isGreenChecked = False

            if not self.showMask:
                rot = np.rot90(self.volArray)
                flip = np.flipud(rot)
                self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            else:
                self.showMaskHelper()

    # hides/shows the probe trajectory
    def toggleProbe(self):
        view = self.image.getView()

        if self.showProbe:
            view.removeItem(self.plItem)
            self.showProbe = False
        else:
            view.addItem(self.plItem)
            self.showProbe = True

    # calculates the euclidean distance between 2 consecutive points annotated on the probe
    def calcDistance(self):
        sorted_dist = sorted(self.distPoints, key=lambda x: x[1])
        for i in range(len(sorted_dist) - 1):
            popup = QMessageBox()
            popup.setText('{} Distance in mm {:.2f}'.format(self.probeDropDown.currentText(), np.linalg.norm(sorted_dist[i] - sorted_dist[i + 1]) / 100))
            popup.exec_()

    # resets the plots
    def resetPlot(self):
        self.plots['unit_density'].resetPlot()
        self.plots[self.metrics.currentText().lower()].resetPlot()

        view = self.image.getView()
        for item in self.lineItems:
            view.removeItem(item)

        for item in self.distItems:
            view.removeItem(item)

        if hasattr(self, 'plItem') and self.prevProbe != self.probeDropDown.currentText():
            view.removeItem(self.plItem)
        
        self.distPoints.clear()
        self.lineItems.clear()
        self.pointsAdded.clear()
        self.oldChannels.clear()

    # hides/shows the probe
    def toggleProbe(self):
        #print(self.show)
        view = self.image.getView()
        if self.show:
            view.removeItem(self.plItem)
            self.show = False
        else:
            view.addItem(self.plItem)
            self.show = True

    # helper function to display the sliced volume that the probe goes through
    def getColorVolume(self, rgb_levels=DEFAULT_COLOR_VALUES):
        level_adjusted_arrays = []
        for colori, int_level in zip(['red', 'green', 'blue'], rgb_levels):
            colarray = np.clip(self.int_arrays[colori], a_min=int_level[0], a_max=int_level[1]) - int_level[0]
            colarray = (colarray * 255. / (int_level[1] - int_level[0])).astype('uint8')
            level_adjusted_arrays.append(colarray)
        return np.stack(level_adjusted_arrays, axis=-1)

    # warps the channels to the ccf
    def warpChannels(self):
        self.channelsOriginal = self.unitPlot.channelsOriginal
        self.channels = self.unitPlot.channels
        values = list(self.acrnm_map.values())
        keys = list(self.acrnm_map.keys())

        # build dataframe for warping
        channel_dict = {'AP': [], 'DV': [], 'ML': [], 'probe_name': []}

        for i in range(len(self.channels)):
            channel = self.channels[i]
            y_coord = int(np.round(channel[1] + 80))

            if y_coord in self.coords:
                coord = self.coords[y_coord] # get the 3d coordinate at that point on the probe track
                channel_dict['AP'].append(coord[0])
                channel_dict['DV'].append(coord[1])
                channel_dict['ML'].append(coord[2])
                

        probe_name = self.probeDropDown.currentText()
        channels = list(range(383, -1, -1))
        probe = [probe_name for i in range(len(channel_dict['AP']))]
        channel_dict['probe_name'] = probe
        df_channel = pd.DataFrame(channel_dict)
        
        # warp using deformation field and reference
        df_final = warp_channels(df_channel, self.field, self.reference, self.probeDropDown.currentText(), channels)
        
        # get areas from ccf
        channel_areas = []

        for index, row in df_final.iterrows():
            struct = self.anno[row.AP, row.DV, row.ML]
            if struct in values:
                ind = values.index(struct)
                key = keys[ind]

                if not key[0].islower():
                    channel_areas.append(key)
                else:
                    channel_areas.append('N/A')
            else:
                channel_areas.append('N/A')

        df_final['channel_areas'] = channel_areas
        df_final.to_csv(os.path.join(self.storageDirectory, '{}_channels_{}_warped.csv'.format(probe_name.replace(' ', '_'), mouse_id)), index=False)

        # group by the region area and plot average for that region
        grouped = df_final.groupby('channel_areas').mean()
        
        for index, row in grouped.iterrows():
            if index != 'N/A':
                plt.text(row.ML, row.DV, index, color='white')
        
        plt.imshow(sitk.GetArrayFromImage(self.reference).T[int(grouped.AP.mean()),0 :, :])
        plt.show()
        

    # function called when the drop down for the metric changes
    # updates the plots to get the metrics file for the corresponding probe and day
    def metricChanged(self):
        print('Metric changed')
        probe = self.probeDropDown.currentText()
        metric = self.metrics.currentText().lower()
        view = self.image.getView()

        probe_let_num = probe[probe.index(' ')+1:]

        key = self.days[int(probe_let_num[1]) - 1]
        paths = self.waveMetricsPath[key]

        self.path = ''
        for p in paths:
            if 'probe' + probe_let_num[0] in p:
                self.path = p
                break

        if self.path != '':
            self.plots[metric].updateMetrics(self.path)
        #self.plots[metric].generateMetrics(metric.lower())

        if hasattr(self, 'linepts'):
            self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues, keep_y=True, old_y=self.plots[self.oldMetric].channels, points_added=self.pointsAdded)

        if hasattr(self, 'oldMetric'):
            view.removeItem(self.plots[self.oldMetric].channelsPlot)
            view.addItem(self.plots[metric].channelsPlot)

        self.oldMetric = metric.lower()

    # displays the region of interest surrounding the probe track
    def displayRegion(self):
        probe = self.probeDropDown.currentText()
        metric = self.metrics.currentText().lower()
        self.path = ''
        view = self.image.getView()

        if probe in self.alignments: # already existing alignment has been done, so use that to update display
            self.alignments[self.prevProbe] = [self.metric, self.plots['unit_density'].channels, self.plots[metric].channels, self.lineItems.copy(), 
                                              self.pointsAdded.copy()] 
            self.resetPlot()
            self.updateDisplay(probe, restore=True)

            # get metrics path for this probe and day
            probe_let_num = probe[probe.index(' ')+1:]

            key = self.days[int(probe_let_num[1]) - 1]
            paths = self.waveMetricsPath[key]

            for p in paths:
                if 'probe' + probe_let_num[0] in p:
                    self.path = p
                    break

            self.plots['unit_density'].updateMetrics(self.path)
            #self.plots['unit_density'].updateDisplay(probe, self.linepts, self.intensityValues) # update display since no existing alignment has been done so far
            self.plots[metric].updateMetrics(self.path)
            self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues, keep_y=True, old_y=self.alignments[probe][1], points_added=self.alignments[probe][-1])

            self.prevProbe = probe
            self.oldMetric = metric.lower()
        else:
            # get metrics path for this probe and day
            probe_let_num = probe[probe.index(' ')+1:]

            key = self.days[int(probe_let_num[1]) - 1]
            paths = self.waveMetricsPath[key]
            for p in paths:
                if 'probe' + probe_let_num[0] in p:
                    self.path = p
                    break

            if self.prevProbe == '' or self.prevProbe != probe: # new alignment
                if self.prevProbe != '' and self.prevProbe != probe:
                    if self.path != '':
                        # add old alignment to dictionary
                        self.alignments[self.prevProbe] = [metric, self.plots['unit_density'].channels, self.plots[metric].channels, self.lineItems.copy(), 
                                                           self.pointsAdded.copy()]
                    else:
                        view.removeItem(self.plots['unit_density'].channelsPlot)
                        view.removeItem(self.plots[metric].channelsPlot)

                    self.resetPlot()
                    self.updateDisplay(probe)
                else: # initial display when nothing has been done
                    self.updateDisplay(probe)

                if self.path != '':
                    self.plots['unit_density'].updateMetrics(self.path)
                    self.plots['unit_density'].updateDisplay(probe, self.linepts, self.intensityValues) # update display since no existing alignment has been done so far

                    self.plots[metric].updateMetrics(self.path)
                    self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues)

                    self.prevProbe = probe
                    self.oldMetric = metric.lower()
                else:
                    popup = QMessageBox()
                    popup.setText('Couldn\'t find metrics.csv for {}'.format(probe))
                    popup.exec_()

    # helper function to remove the alignment line
    def removeLineHelper(self, y_coord, ind):
        view = self.image.getView()
        item = self.lineItems[ind]
        view.removeItem(item)
        self.lineItems.remove(item)
        self.pointsAdded.remove(y_coord)

        self.plots['unit_density'].channelsPlot.setData(pos=np.array(self.plots['unit_density'].channels, dtype=float), adj=np.array(self.plots['unit_density'].adj, dtype=int))
        self.plots[self.metrics.currentText().lower()].channelsPlot.setData(pos=np.array(self.plots[self.metrics.currentText().lower()].channels, dtype=float), 
                                                                            adj=np.array(self.plots[self.metrics.currentText().lower()].adj, dtype=int))
        self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
        self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)


    # remove the line clicked on from the refinement display
    def removeLine(self, plots, points):
        y_coord = points[0].pos()[1]
        ind = self.pointsAdded.index(y_coord)
        self.removeLineHelper(y_coord, ind)

    # updates the other plots to align it with the plot that was clicked on (the unit density plot in this case)
    def updateAfterClick(self, other_plot, line_point, channel, pointsAdded):
        #print('Other plot clicked', other_plot.measurement)
        other_plot.channelsPlot.scatterPoint = other_plot.channels[channel]
        other_plot.oldChannels.append(other_plot.channels)
        flag = True

        if line_point[1] != other_plot.channelsPlot.scatterPoint[1]:
            if line_point[1] < other_plot.channelsPlot.scatterPoint[1]:
                #print('lower')
                diff = abs(line_point[1] - other_plot.channelsPlot.scatterPoint[1])
                newPoints = []

                if len(pointsAdded) == 0:
                    for p in other_plot.channels:
                        newPoints.append([p[0], p[1] - diff])
                else:
                    srt = sorted(pointsAdded, reverse=True)
                    greater = [t for t in pointsAdded if t > line_point[1]]
                    less = [t for t in pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 exisiting, don't shift, just linearly interpolate
                        flag = other_plot.linspaceAbove(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new line is below exiting, so shift
                        for p in other_plot.channels:
                            if p[1] <= other_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else: # new alignment is above, so interpolate and shift
                        flag = other_plot.linspaceAbove(srt, channel, diff, flag)

                        for i in range(len(other_plot.channels)):
                            p = other_plot.channels[i]
                            if p[1] > other_plot.channelsPlot.scatterPoint[1]:
                                other_plot.channels[i] = [p[0], p[1] - diff]
            elif line_point[1] > other_plot.channelsPlot.scatterPoint[1]:
                diff = line_point[1] - other_plot.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(pointsAdded) == 0:
                    for p in other_plot.channels:
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(pointsAdded)
                    greater = [t for t in pointsAdded if t > line_point[1]]
                    less = [t for t in pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0:
                        flag = other_plot.linspaceBelow(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0:
                        flag = other_plot.linspaceBelow(srt, channel, diff, flag)

                        for i in range(len(other_plot.channels)):
                            p = other_plot.channels[i]
                            if p[1] < other_plot.channelsPlot.scatterPoint[1]:
                                other_plot.channels[i] = [p[0], p[1] + diff]
                    else:
                        for p in other_plot.channels:
                            if p[1] >= other_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
            
            if flag:
                other_plot.channels = newPoints
            #if len(self.pointsAdded) == 1:
            other_plot.channelsPlot.setData(pos=np.array(other_plot.channels, dtype=float), adj=np.array(other_plot.adj, dtype=int))

        #other_plot.linearSpacePoints(pointsAdded)

    # pefroms alignment on the unit density plot
    def onClickProbeHelper(self, click_plot, line_point, is_unit=True):
        channel = click_plot.channels.index([click_plot.channelsPlot.scatterPoint[0], click_plot.channelsPlot.scatterPoint[1]]) # get index of point clicked on in unit density
        click_plot.oldChannels.append(click_plot.channels)
        flag = True

        if line_point[1] != click_plot.channelsPlot.scatterPoint[1]:
            if line_point[1] < click_plot.channelsPlot.scatterPoint[1]: # alignment is above currently clicked channel
                #print('lower')
                diff = abs(line_point[1] - click_plot.channelsPlot.scatterPoint[1])
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in click_plot.channels:
                        newPoints.append([p[0], p[1] - diff])
                else:
                    srt = sorted(self.pointsAdded, reverse=True)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 exisiting, don't shift, just linearly interpolate
                        flag = click_plot.linspaceAbove(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new line is below existing, so shift
                        for p in click_plot.channels:
                            if p[1] <= click_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else: # new line is above so linearly interpolate and shift
                        flag = click_plot.linspaceAbove(srt, channel, diff, flag)

                        for i in range(len(click_plot.channels)):
                            p = click_plot.channels[i]
                            if p[1] > click_plot.channelsPlot.scatterPoint[1]:
                                click_plot.channels[i] = [p[0], p[1] - diff]
            elif line_point[1] > click_plot.channelsPlot.scatterPoint[1]: # alignment is below currently clicked channel
                diff = line_point[1] - click_plot.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in click_plot.channels:
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(self.pointsAdded)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 existing, don't shift, just linearly interpolate
                        flag = click_plot.linspaceBelow(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new alignment is above, so linearly interpolate and shift
                        flag = click_plot.linspaceBelow(srt, channel, diff, flag)

                        for i in range(len(click_plot.channels)):
                            p = click_plot.channels[i]
                            if p[1] < click_plot.channelsPlot.scatterPoint[1]:
                                click_plot.channels[i] = [p[0], p[1] + diff]
                    else: # new alignment is below, so just shift
                        for p in click_plot.channels:
                            if p[1] >= click_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
                        
            if flag:
                click_plot.channels = newPoints
            #if len(self.pointsAdded) == 1:
            click_plot.channelsPlot.setData(pos=np.array(click_plot.channels, dtype=float), adj=np.array(click_plot.adj, dtype=int))
        
            pts = [[t, line_point[1]] for t in range(int(click_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)
        else: # alignment is same as y coord of the channel clicked on
            pts = [[t, line_point[1]] for t in range(int(click_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            diff = 0
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)

        # update other plots to match alignment of unit density
        if is_unit:
            self.updateAfterClick(self.plots[self.metrics.currentText().lower()], line_point, channel, self.pointsAdded)
        else:
            self.updateAfterClick(self.plots['unit_density'], line_point, channel, self.pointsAdded)

        self.pointsAdded.append(line_point[1])
        click_plot.linearSpacePoints(self.pointsAdded)
        
        if is_unit:
            self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)
        else:
            self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
        
        view = self.image.getView()
        view.addItem(h_line)
        #h_line.setClickable(True)
        h_line.sigClicked.connect(self.removeLine)

    def onclickProbe(self, plot, points):
        if self.plots['unit_density'].channelsPlot.clicked:
            line_point = points[0].pos()
            self.onClickProbeHelper(self.plots['unit_density'], line_point, is_unit=True)
        else:
            line_point = points[0].pos()
            self.onClickProbeHelper(self.plots[self.metrics.currentText().lower()], line_point, is_unit=False)

        self.plots['unit_density'].channelsPlot.clicked = False
        self.plots[self.metrics.currentText().lower()].channelsPlot.clicked = False
    
    # left arrow key
    # removes the most recent alignment
    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Left:
            popped = self.lineItems.pop(-1) # last alignment done
            view = self.image.getView()
            view.removeItem(popped)
            self.pointsAdded.pop(-1)
            self.plots['unit_density'].linearSpacePoints(self.pointsAdded)
            self.plots[self.metrics.currentText().lower()].linearSpacePoints(self.pointsAdded)

    # toggles the mask 
    def toggleCCFRegions(self):
        if self.showMask:
            rot = np.rot90(self.volArray)
            flip = np.flipud(rot)
            self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
            self.showMask = False
        else:
            self.showMaskHelper()
            self.showMask = True

    # displays the region along the probe track
    # probe: string, the probe to be displayed from the drop down
    def updateDisplay(self, probe, restore=False):
        self.showMask = False
        self.showProbe = True
        #print('Probe', probe)
        x = self.probeAnnotations[self.probeAnnotations.probe_name == probe].ML 
        y = self.probeAnnotations[self.probeAnnotations.probe_name == probe].DV 
        
        z = self.probeAnnotations[self.probeAnnotations.probe_name == probe].AP

        # get trajectory
        if len(z) > 0: # 
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
            
        print(linepts.shape[0]) 
        intensity_values_red = np.zeros((linepts.shape[0],160)) # dummy array to be used later on for class fields below
        self.coords = {}
        self.linepts = linepts
        self.intensityValues = intensity_values_red

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2]) # dictionary to store 3d coord at probe y coord

        p = probe.replace(' ', '_')
        self.volArray = np.array(Image.open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_slice.png'.format(p)))) # read slice 
        self.mask = Image.open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_mask.png'.format(p))) # read mask
        self.blended = np.array(Image.open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_overlay.png'.format(p)))) # read overlay
        with open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_labels.pickle'.format(p)), 'rb') as handle:
            self.labels_pos = pickle.load(handle)

        # modify based on if red and green channel toggle have been checked or not
        if self.isRedChecked:
            self.redOld = self.volArray[:, :, 0].copy()
            self.volArray[:, :, 0] = 0

        if self.isGreenChecked:
            self.greenOld = self.volArray[:, :, 1].copy()
            self.volArray[:, :, 1] = 0

        
        #self.volArray = np.array(vol_mask)
        #self.volArray *= 2
        rot = np.rot90(self.volArray)
        flip = np.flipud(rot)
        rot_mask = np.rot90(self.mask)
        flip_mask = np.flipud(rot_mask)

        self.imageMask.setImage(flip_mask, levels=(0, 255), autoRange=False)
        self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
        view = self.image.getView()

        self.points = [[80, t] for t in range(j)]
        self.plItem = pg.ScatterPlotItem(pos=self.points, pen=QtGui.QPen(QColor('red')), brush=QtGui.QBrush(QColor('red')))
        #self.plItem.setClickable(True)
        self.plItem.sigClicked.connect(self.onclickProbe)

        if not restore:
            if self.path != '':
                view.addItem(self.plots['unit_density'].channelsPlot)
                view.addItem(self.plots[self.metrics.currentText().lower()].channelsPlot)

            view.addItem(self.plItem)
        else: # read from dictionary storing saved plots for the probe
            view.addItem(self.plItem)
            plot_items = self.alignments[probe]

            self.plots['unit_density'].channelsPlot.setData(pos=np.array(plot_items[1], dtype=float), adj=np.array(self.plots['unit_density'].adj, dtype=int))
            self.plots['unit_density'].channels = plot_items[1].copy()
            self.plots[self.metrics.currentText().lower()].channelsPlot.setData(pos=np.array(plot_items[2], dtype=float), adj=np.array(self.plots['unit_density'].adj, dtype=int))
            self.plots[self.metrics.currentText().lower()].channels = plot_items[2].copy()


            view.addItem(self.plots['unit_density'].channelsPlot) # add unit plot
            view.addItem(self.plots[self.metrics.currentText().lower()].channelsPlot) # add metric plot

            self.lineItems = plot_items[3].copy() # add line items
            for item in self.lineItems:
                view.addItem(item)

            self.pointsAdded = plot_items[4].copy() # restore alignment

        #view.addItem(self.textItem)
        
if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID

    app = QApplication(sys.argv)
    v = VolumeAlignment(mouse_id)
    sys.exit(app.exec_())