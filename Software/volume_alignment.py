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

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)

SCALING_FACTOR = 1.5
DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]

class Graph(pg.GraphItem):
    def __init__(self):
        self.dragPoint = None
        self.dragOffset = None
        self.textItems = []
        self.pointPos = {}
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

        points[0].setBrush(QtGui.QBrush(QColor('green')))
        points[0].setPen(QtGui.QPen(QColor('green')))
        self.lastPoint = points[0]
        self.scatterPoint = points[0].pos()

class VolumeAlignment(QWidget):
    def __init__(self, mouse_id):
        super().__init__()
        self.title = 'Volume Alignment'
        self.mouseID = mouse_id
        self.unitDense = pd.read_csv(pathlib.Path('//allen/programs/mindscope/workgroups/dynamicrouting/PilotEphys/Task 2 pilot/2022-08-15_11-22-28_626791/Record Node 108/experiment1/recording1/continuous/Neuropix-PXI-102.ProbeA-AP/waveform_metrics.csv'))

        self.initUI()
    
    def initUI(self):
        self.setWindowTitle(self.title)
        self.width = int(4000) #default = coronal view
        self.height = int(4000)

        # create default image
        im8 = np.ones((self.height,self.width),dtype='uint8')*255
        self.image = pg.image(im8)
        self.image.ui.histogram.hide()
        self.image.ui.roiBtn.hide()
        self.image.ui.menuBtn.hide()
        self.image.setObjectName('image')
        self.image.setImage(im8.transpose())
        self.pointsAdded = []
        self.channelCoords = {}
        self.lineItems = []
        self.oldChannels = [] # stack for undoing lines

        #self.unit_dense = np.random.randint(384, size=384)
        self.frequencyCounts = self.unitDense['peak_channel'].value_counts().sort_index().reset_index().to_numpy()
        #x = np.linspace(-10, 100, num=384)
        self.channelsRecorded = self.frequencyCounts[:, 0].tolist()
        self.channelsOriginal = []

        for i in range(384):
            if i in self.channelsRecorded:
                ind = self.channelsRecorded.index(i)
                self.channelsOriginal.append([i, self.frequencyCounts[ind, 1]])
            else:
                self.channelsOriginal.append([i, 0])

        self.channelsOriginal = [[self.channelsOriginal[i][1], self.channelsOriginal[len(self.channelsOriginal) - i - 1][0]] for i in range(len(self.channelsOriginal))]
        #self.densityChannels = [[50, int(i * 3.84)] for i in range(384)]
        self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]

        # metric plot
        self.channelsPlot = Graph()
        #self.denseChannelsPlot = Graph()

        #self.denseChannelsPlot.setData(pos=np.array(self.densityChannels, dtype=int), adj=np.array(self.adj, dtype=int))
        self.channelsPlot.setData(pos=np.array(self.channelsOriginal, dtype=float), adj=np.array(self.adj, dtype=int))

        #view = self.image.getView()
        #view.addItem(self.channelsPlot)
        #view.addItem(self.denseChannelsPlot)

        # important directories
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_file)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')

        self.reference = sitk.ReadImage( self.reference_file )
        self.field = sitk.ReadImage( self.field_file )

        self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))
        self.volumeImage = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_green.mhd'))).T
        #self.getColorVolume()

        # main layout with region of interest
        self.mainLayout = QVBoxLayout()
        self.mainLayout.addWidget(self.image)

        self.probes = self.probeAnnotations['probe_name'].unique()
        self.probeDropDown = QComboBox()
        for probe in self.probes:
            self.probeDropDown.addItem(probe)

        self.metrics = QComboBox()
        self.metrics.addItem('Metrics')
        self.metrics.addItem('Unit Density')

        self.viewButton = QPushButton('View Probe Region with Selected Metric')
        self.viewButton.clicked.connect(self.displayRegion)

        self.resetPlotButton = QPushButton('Reset Metric Plot')
        self.resetPlotButton.clicked.connect(self.resetPlot)

        self.warpButton = QPushButton('Warp to CCF')
        self.warpButton.clicked.connect(self.warpChannels)

        self.probeViewLayout = QHBoxLayout()
        self.probeViewLayout.addWidget(self.probeDropDown)
        self.probeViewLayout.addWidget(self.metrics)
        self.probeViewLayout.addWidget(self.viewButton)
        self.probeViewLayout.addWidget(self.resetPlotButton)

        self.probeViewLayout.addWidget(self.warpButton)
        self.probeViewLayout.setAlignment(QtCore.Qt.AlignTop)
        self.mainLayout.addLayout(self.probeViewLayout)

        self.setLayout(self.mainLayout)
        self.showMaximized()


    # when a new probe is displayed
    def resetPlot(self):
        self.channels = self.channelsOriginal
        self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

        view = self.image.getView()
        for item in self.lineItems:
            view.removeItem(item)

        self.lineItems.clear()
        self.pointsAdded.clear()
        self.oldChannels.clear()

    def getColorVolume(self, intensity_values, rgb_levels=DEFAULT_COLOR_VALUES):
       colarray = np.clip(intensity_values, a_min=rgb_levels[0][0], a_max=rgb_levels[0][1]) - rgb_levels[0][0]
       colarray = (colarray * 255. / (rgb_levels[0][1] - rgb_levels[0][0])).astype('uint8')
       
       return colarray

    # warps the channels to the ccf
    def warpChannels(self):
        channel_dict = {'AP': [], 'DV': [], 'ML': [], 'probe_name': []}
        channels = [p[1] for p in self.channelsOriginal]

        for i in range(len(self.channels)):
            channel = self.channels[i]
            y_coord = int(np.round(channel[1]))

            if y_coord not in self.coords:
                coord = (0, 0, 0)
            else:
                coord = self.coords[y_coord] # get the 3d coordinate at that point on the probe track

            print(coord)
            channel_dict['AP'].append(coord[0])
            channel_dict['DV'].append(coord[1])
            channel_dict['ML'].append(coord[2] * 0.875)

        probe_name = self.probeDropDown.currentText()
        probe = [probe_name for i in range(len(channel_dict['AP']))]
        channel_dict['probe_name'] = probe
        df_channel = pd.DataFrame(channel_dict)
        
        warp_channels(self.storageDirectory, df_channel, self.field, self.reference, self.probeDropDown.currentText(), self.mouseID, channels)

    # displays the region of interest surrounding the probe track
    def displayRegion(self):
        probe = self.probeDropDown.currentText()
        view = self.image.getView()
        self.resetPlot()
        self.updateDisplay(probe)
        view.addItem(self.channelsPlot)

    # helper function for shfiting points when new alignment is added
    def closest(self, lst, K):
        if len(lst) > 0:
            return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

    # remove the line clicked on from the refinement display
    def removeLine(self, plots, points):
        y_coord = points[0].pos()[1]
        ind = self.pointsAdded.index(y_coord)

        view = self.image.getView()
        item = self.lineItems[ind]
        view.removeItem(item)
        self.lineItems.remove(item)
        self.pointsAdded.remove(y_coord)

        if len(self.pointsAdded) == 0:
            self.channels = self.ogPoints

        self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        self.linearSpacePoints()

    def replaceValues(self, lp, points_between):
        #print(points_between)
        for i in range(len(lp)):
            ind = self.channels.index(points_between[len(lp) - i - 1])
            #print(self.channels[ind])
            self.channels[ind][1] = lp[i]
            #print(self.channels[ind])

    def linearSpacePoints(self):
        if len(self.pointsAdded) > 1:
            sorted_linepts = sorted(self.pointsAdded)

            for i in range(len(sorted_linepts) - 1): # each time new line is added
                anchor_top = sorted_linepts[i]
                anchor_bottom = sorted_linepts[i + 1]

                #print(anchor_top, anchor_bottom)

                points_between = [p for p in self.channels if p[1] > anchor_top and p[1] < anchor_bottom]
                #print(len(points_between))
                lp = np.linspace(anchor_top + 1, anchor_bottom - 1, num=len(points_between)).tolist()
                #print(len(lp))
                #print(lp)
                #lp = [[p[0], int(p[1])] for p in lp]
                self.replaceValues(lp, points_between)
                #print(len(self.channels))

            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

    def onclickProbe(self, plot, points):
        print(len(self.channels))
        line_point = points[0].pos()
        #print(self.channels)
        channel = self.channels.index([self.channelsPlot.scatterPoint[0], self.channelsPlot.scatterPoint[1]])

        self.oldChannels.append(self.channels)
        flag = True

        if line_point[1] != self.channelsPlot.scatterPoint[1]:
            if line_point[1] < self.channelsPlot.scatterPoint[1]:
                #print('lower')
                diff = abs(line_point[1] - self.channelsPlot.scatterPoint[1])
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in self.channels:
                        """
                        if p[1] <= self.channelsPlot.scatterPoint[1]:
                            newPoints.append([p[0], p[1] - diff])
                        else:
                            newPoints.append([p[0], p[1]])
                        """
                        newPoints.append([p[0], p[1] - diff])
                else:
                    srt = sorted(self.pointsAdded, reverse=True)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 exisiting, don't shift, just linearly interpolate
                        for point in srt:
                            if self.channelsPlot.scatterPoint[1] > point:
                                points_between = [p for p in self.channels if p[1] < self.channelsPlot.scatterPoint[1] and p[1] > point]
                                self.channels[channel][1] -= diff
                                lp = np.linspace(point + 1, self.channels[channel][1] - 1, len(points_between))
                                self.replaceValues(lp, points_between)
                                flag = False
                                break
                    elif len(greater) > 0 and len(less) == 0:
                        for p in self.channels:
                            if p[1] <= self.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else:
                        for point in srt:
                            if self.channelsPlot.scatterPoint[1] > point:
                                points_between = [p for p in self.channels if p[1] < self.channelsPlot.scatterPoint[1] and p[1] > point]
                                self.channels[channel][1] -= diff
                                lp = np.linspace(point + 1, self.channels[channel][1] - 1, len(points_between))
                                self.replaceValues(lp, points_between)
                                flag = False
                                break

                """
                for p in self.channels:
                    if len(self.pointsAdded) == 0:
                        if p[1] <= self.channelsPlot.scatterPoint[1]:
                            newPoints.append([p[0], p[1] - diff])
                        else:
                            newPoints.append([p[0], p[1]])
                    else:
                        diff_closest = p[1] - self.closest(self.pointsAdded, p[1])

                        if diff_closest > 0:
                            if p[1] <= self.channelsPlot.scatterPoint[1] and p[1] > self.closest(self.pointsAdded, p[1]) and p[1] not in self.pointsAdded:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                        else:
                            if p[1] <= self.channelsPlot.scatterPoint[1] and p[1] < self.closest(self.pointsAdded, p[1]) and p[1] not in self.pointsAdded:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                """
            elif line_point[1] > self.channelsPlot.scatterPoint[1]:
                diff = line_point[1] - self.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in self.channels:
                        """
                        if p[1] >= self.channelsPlot.scatterPoint[1]:
                            newPoints.append([p[0], p[1] + diff])
                        else:
                            newPoints.append([p[0], p[1]])
                        """
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(self.pointsAdded)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0:
                        for point in srt:
                            if self.channelsPlot.scatterPoint[1] < point:
                                points_between = [p for p in self.channels if p[1] > self.channelsPlot.scatterPoint[1] and p[1] < point]
                                self.channels[channel][1] += diff
                                lp = np.linspace(self.channels[channel][1] + 1, point - 1, len(points_between))
                                #print(lp)
                                self.replaceValues(lp, points_between)
                                flag = False
                                break
                    elif len(greater) > 0 and len(less) == 0:
                        for point in srt:
                            if self.channelsPlot.scatterPoint[1] < point:
                                points_between = [p for p in self.channels if p[1] > self.channelsPlot.scatterPoint[1] and p[1] < point]
                                self.channels[channel][1] += diff
                                lp = np.linspace(self.channels[channel][1] + 1, point - 1, len(points_between))
                                #print(lp)
                                self.replaceValues(lp, points_between)
                                flag = False
                                break
                    else:
                        for p in self.channels:
                            if p[1] >= self.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])

                        """
                        for p in self.channels:
                            
                            if p[1] + diff >= line_point[1] and p[1] not in self.pointsAdded:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
                            
                            if diff_closest > 0:
                                if p[1] + diff  >= line_point[1] and p[1] + diff > self.closest(self.pointsAdded, p[1]) and p[1] not in self.pointsAdded:
                                    newPoints.append([p[0], p[1] + diff])
                                else:
                                    newPoints.append([p[0], p[1]])
                            else:
                                if p[1] <= self.channelsPlot.scatterPoint[1] and p[1] + diff < self.closest(self.pointsAdded, p[1]) and p[1] not in self.pointsAdded:
                                    newPoints.append([p[0], p[1] + diff])
                                else:
                                    newPoints.append([p[0], p[1]])
                            
                    
                    if len(self.pointsAdded) == 0:
                        if p[1] >= self.channelsPlot.scatterPoint[1]:
                            newPoints.append([p[0], p[1] + diff])
                        else:
                            newPoints.append([p[0], p[1]])
                    else:
                        diff_closest = p[1] - self.closest(self.pointsAdded, p[1])

                        if diff_closest > 0:
                            if p[1] >= self.channelsPlot.scatterPoint[1] and p[1] > self.closest(self.pointsAdded, p[1]) and p[1] not in self.pointsAdded:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
                        else:
                            if p[1] >= self.channelsPlot.scatterPoint[1] and p[1] < self.closest(self.pointsAdded, p[1]) and p[1] not in self.pointsAdded:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    """
            
            if flag:
                self.channels = newPoints
            #if len(self.pointsAdded) == 1:
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        
            pts = [[t, line_point[1]] for t in range(int(self.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            #h_line = pg.PlotCurveItem([self.channelsPlot.scatterPoint[0], line_point[0]], [line_point[1], line_point[1]], pen=QtGui.QPen(QColor('yellow')),
                                          #brush=QtGui.QBrush(QColor('yellow')))
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)
        else:
            pts = [[t, line_point[1]] for t in range(int(self.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            diff = 0
            #h_line = pg.PlotCurveItem([self.channelsPlot.scatterPoint[0], line_point[0]], [self.channelsPlot.scatterPoint[1], line_point[1]], pen=QtGui.QPen(QColor('yellow')),
                                        #brush=QtGui.QBrush(QColor('yellow')))
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)

        #self.pointsAdded.append(line_point[1])

        if line_point[1] in self.coords:
            self.channelCoords[channel] = self.coords[line_point[1]]

        self.pointsAdded.append(line_point[1])
   
        self.linearSpacePoints()
        #print(self.channelCoords)
        view = self.image.getView()
        view.addItem(h_line)
        #h_line.setClickable(True)
        h_line.sigClicked.connect(self.removeLine)

    def getAffineMatrix(self, direction='tvr'):
        storage_directory = get_tc_info(self.mouseID)

        with open('/{}/local_alignment/localalignment_transform_input.xml'.format(storage_directory)) as f:
            transform_dict = xmltodict.parse(f.read())
            align3d = transform_dict['alignment']['alignment3ds']['alignment3d']
            print(align3d)

            self.affine = np.zeros((4, 3), dtype=float)
            i = 0
            j = 0

            for key in align3d:
                if direction in key:
                    self.affine[i][j] = float(align3d[key]['#text'])

                    if j == self.affine.shape[1] - 1: # move to next row in matrix
                        j = 0
                        i += 1
                    else:
                        j += 1

                if i == self.affine.shape[0]: # finished
                    break

    # displays the region along the probe track
    # probe: the probe to be displayed from the drop down
    def updateDisplay(self, probe):
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

            linepts = vv[0] * np.mgrid[-530:600:0.7][:,np.newaxis]
            linepts += datamean
            
            if linepts[-1,1] - linepts[0,1] < 0:
                linepts = np.flipud(linepts)
            
        intensity_values = np.zeros((linepts.shape[0],160))
        self.coords = {}
        self.linepts = linepts

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2])

            for k in range(-80,80):
                try:
                    intensity_values[j,k+80] = (self.volumeImage[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                except IndexError:
                    pass
        
        # display image
        #intensity_values *= 2.
        #self.trigScale = intensity_values.shape[0] / intensity_values.shape[1] 
        
        self.trigScale = (np.linalg.norm(linepts[linepts.shape[0] - 1, :]) - np.linalg.norm(linepts[0, :])) / intensity_values.shape[1]
        print('Trig Scale', self.trigScale)
        self.getAffineMatrix()

        self.dz = linepts[linepts.shape[0] - 1, 0] - linepts[0, 0]
        self.dy = linepts[linepts.shape[0] - 1, 1] - linepts[0, 1]
        self.dx = linepts[linepts.shape[0] - 1, 2] - linepts[0, 2]
        self.vector = np.array([self.dz, self.dy, self.dx])
        self.vectorPrime = self.affine.dot(self.vector)
        self.affScale = np.linalg.norm(self.vectorPrime[0:3]) / np.linalg.norm(self.vector)
        print('Aff Scale', self.affScale)

        self.scale = 3840 / (self.trigScale * 10 *  (1 / self.affScale))
        print('Scale', self.scale)
        self.scale /= 384

        print('Final Scale', self.scale)
        newPoints = []
        for p in self.channelsOriginal:
            if self.scale >= 1:
                newPoints.append([p[0], ((p[1] / self.scale) + p[1])])
            else:
                newPoints.append([p[0], ((p[1] * self.scale) + p[1])])
        
        self.ogPoints = newPoints
        self.channels = newPoints
        self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        
        self.volArray = self.getColorVolume(intensity_values)
        print(self.volArray.shape)
        
        self.volArray *= 2
        rot = np.rot90(self.volArray)
        flip = np.flipud(rot)
      
        self.image.setImage(flip[:, 500:j], levels=(0, 255), autoRange=False)

        view = self.image.getView()
        self.points = [[80, t] for t in range(-500, j)]
        self.plItem = pg.ScatterPlotItem(pos=self.points, pen=QtGui.QPen(QColor('red')), brush=QtGui.QBrush(QColor('red')))
        #self.plItem.setClickable(True)
        self.plItem.sigClicked.connect(self.onclickProbe)
        view.addItem(self.plItem)

if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID

    app = QApplication(sys.argv)
    w = VolumeAlignment(mouse_id)
    sys.exit(app.exec_())