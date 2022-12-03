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

class PlotDisplayItem():
    # create default image
    def __init__(self, measurement, waveform, metrics, volume_image, probe_annotations, mouse_id, metrics_list):
        self.width = int(4000) #default = coronal view
        self.height = int(4000)
        self.remove = True
        self.show = True

        self.probeAnnotations = probe_annotations
        self.mouseID = mouse_id
        """
        im8 = np.ones((self.height,self.width),dtype='uint8')*255
        self.image = pg.image(im8)
        self.image.ui.histogram.hide()
        self.image.ui.roiBtn.hide()
        self.image.ui.menuBtn.hide()
        self.image.setObjectName('image')
        #self.image.setImage(im8.transpose())
        """
        self.pointsAdded = []
        self.channelCoords = {}
        self.lineItems = []
        self.oldChannels = [] # stack for undoing lines
        self.waveform = waveform
        self.metrics = metrics
        self.metricsList = metrics_list

        self.waveform_metrics = self.waveform.merge(self.metrics, on='cluster_id')
        self.measurement = measurement
        self.otherPlots = []
        #self.generateMetrics(measurement)

        self.volumeImage = volume_image
        self.channelsPlot = Graph()
        #self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
        self.textItem = pg.TextItem(measurement.upper(), anchor=(1, 1))

        #self.channels = self.channelsOriginal
    
    def updateMetrics(self, path):
        self.waveform_metrics = pd.read_csv(path)
        print('Metrics Path', path)
        self.generateMetrics(self.measurement)

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
           
            self.channelsOriginal = [[self.channelsOriginal[i][1] - 20, i] for i in range(len(self.channelsOriginal))]
        elif self.measurement == 'spread':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=0, rolling_value=150)
        elif self.measurement == 'firing_rate':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=3, rolling_value=150)
        elif self.measurement == 'd_prime':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=1, rolling_value=150)
        elif self.measurement == 'cumulative_drift':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=5, rolling_value=150)
        elif self.measurement == 'velocity_above':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=1/2, rolling_value=150)
        elif self.measurement == 'velocity_below':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=1/2, rolling_value=150)
        elif self.measurement == 'amplitude':
            self.processMetrics()
            self.generateMetricChannels(measurement, shift_value=5, rolling_value=150)

    def generateMetricChannels(self, metric, shift_value, rolling_value):
        print('Metric', metric)
        peak_values = self.averageMetricsChannels['peak_channel'].values.tolist()
        values = self.averageMetricsChannels[metric].to_numpy().tolist()

        if 'velocity' in metric:
            conv = np.ones(2)
        else:
            conv = np.ones(20)
        #peak_values = [int(p) for p in peak_values]
        self.channelsOriginal = []

        for i in range(384):
            if i in peak_values:
                index = peak_values.index(i)
                self.channelsOriginal.append([values[index], i])
                #conv[index] = 1
            else:
                self.channelsOriginal.append([0, i])

        x_val = [p[0] for p in self.channelsOriginal]
        smoothed = np.convolve(x_val, conv, mode='same') / np.sum(conv)
        print(smoothed.shape)
        for i in range(384):
            if shift_value != 0:
                self.channelsOriginal[i] = [(smoothed[i] / shift_value) - rolling_value, self.channelsOriginal[i][1]]
            else:
                self.channelsOriginal[i] = [(smoothed[i]) - rolling_value, self.channelsOriginal[i][1]]

    def processMetrics(self):
        self.averageMetricsChannels = self.waveform_metrics.groupby('peak_channel').mean().reset_index()
        #self.averageMetricsChannels = self.averageMetricsChannels.rolling(rolling_value, win_type='boxcar').mean()
    # when a new probe is displayed
    def resetPlot(self, remove_probe=False):
        self.channels = self.channelsOriginal
        self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

    # helper function for shfiting points when new alignment is added
    def closest(self, lst, K):
        if len(lst) > 0:
            return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

    """
    # function to remove lines from other plots to keep alignment across all plots
    def removeLineOthers(self, other_plot, ind):
        y_coord = other_plot.pointsAdded[ind]

        view = other_plot.image.getView()
        item = other_plot.lineItems[ind]
        view.removeItem(item)
        other_plot.lineItems.remove(item)
        other_plot.pointsAdded.remove(y_coord)

        if len(other_plot.pointsAdded) == 0:
            other_plot.channels = other_plot.ogPoints

        other_plot.channelsPlot.setData(pos=np.array(other_plot.channels, dtype=float), adj=np.array(other_plot.adj, dtype=int))
        other_plot.linearSpacePoints()

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

        for plot in self.otherPlots:
            self.removeLineOthers(plot, ind)
    """
    def replaceValues(self, lp, points_between):
        #print(points_between)
        for i in range(len(lp)):
            ind = self.channels.index(points_between[i])
            #print(self.channels[ind])
            self.channels[ind][1] = lp[i]
            #print(self.channels[ind])

    def linearSpacePoints(self, pointsAdded):
        if len(pointsAdded) > 1:
            sorted_linepts = sorted(pointsAdded)

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

                print('Anchor Top', anchor_top)
                print('Anchor Bottom', anchor_bottom)
                scale = (anchor_bottom - anchor_top) / len(points_between)
                print('Scale', scale)

            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))

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
    """
    # updates the other plots to align it with the plot that was clicked on
    def updateAfterClick(self, other_plot, line_point, channel, pointsAdded):
        print('Other plot clicked', other_plot.measurement)
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
                        flag = other_plot.linspaceAbove(srt, channel, diff, flag, other_plot)
                    elif len(greater) > 0 and len(less) == 0: # new line is below exiting, so shift
                        for p in other_plot.channels:
                            if p[1] <= other_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else:
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

        other_plot.linearSpacePoints(pointsAdded)
        
            pts = [[t, line_point[1]] for t in range(int(other_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            other_plot.lineItems.append(h_line)
        else:
            pts = [[t, line_point[1]] for t in range(int(other_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            diff = 0
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            other_plot.lineItems.append(h_line)

        #self.pointsAdded.append(line_point[1])
        
        if line_point[1] in other_plot.coords:
            other_plot.channelCoords[channel] = other_plot.coords[line_point[1]]

        other_plot.pointsAdded.append(line_point[1])
        #print(other_plot.coords[int(line_point[1])])
        other_plot.linearSpacePoints()
        
        #print(self.channelCoords)
        view = other_plot.image.getView()
        view.addItem(h_line)
        #h_line.setClickable(True)
        h_line.sigClicked.connect(other_plot.removeLine)
        """

    """
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
                        newPoints.append([p[0], p[1] - diff])
                else:
                    srt = sorted(self.pointsAdded, reverse=True)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0: # new alignment is between 2 exisiting, don't shift, just linearly interpolate
                        flag = self.linspaceAbove(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0: # new line is below exiting, so shift
                        for p in self.channels:
                            if p[1] <= self.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else:
                        flag = self.linspaceAbove(srt, channel, diff, flag)

                        for i in range(len(self.channels)):
                            p = self.channels[i]
                            if p[1] > self.channelsPlot.scatterPoint[1]:
                                self.channels[i] = [p[0], p[1] - diff]
            elif line_point[1] > self.channelsPlot.scatterPoint[1]:
                diff = line_point[1] - self.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in self.channels:
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(self.pointsAdded)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0:
                        flag = self.linspaceBelow(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0:
                        flag = self.linspaceBelow(srt, channel, diff, flag)

                        for i in range(len(self.channels)):
                            p = self.channels[i]
                            if p[1] < self.channelsPlot.scatterPoint[1]:
                                self.channels[i] = [p[0], p[1] + diff]
                    else:
                        for p in self.channels:
                            if p[1] >= self.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] + diff])
                            else:
                                newPoints.append([p[0], p[1]])
            
            if flag:
                self.channels = newPoints
            #if len(self.pointsAdded) == 1:
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        
            pts = [[t, line_point[1]] for t in range(int(self.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)
        else:
            pts = [[t, line_point[1]] for t in range(int(self.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            diff = 0
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)

        #self.pointsAdded.append(line_point[1])

        if line_point[1] in self.coords:
            self.channelCoords[channel] = self.coords[line_point[1]]

        self.pointsAdded.append(line_point[1])
        print(self.coords[int(line_point[1])])
        self.linearSpacePoints()
        
        #print(self.channelCoords)
        view = self.image.getView()
        view.addItem(h_line)
        #h_line.setClickable(True)
        h_line.sigClicked.connect(self.removeLine)

        for plot in self.otherPlots:
            print(plot.measurement)
            self.updateAfterClick(plot, line_point, channel)
    """
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

    # displays the region along the probe track
    # probe: the probe to be displayed from the drop down
    def updateDisplay(self, probe, linepts, intensity_values, keep_y=False, old_y=None, points_added=None):
        """
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
        """
        # display image
        #intensity_values *= 2.
        #self.trigScale = intensity_values.shape[0] / intensity_values.shape[1] 
        
        #print(intensity_values.shape)
        
        if not keep_y:
            self.trigScale = np.linalg.norm(linepts[-1] - linepts[0]) / intensity_values.shape[0]
            #self.trigScale = 1
            print('Trig Scale', self.trigScale)
            self.getAffineMatrix()
            print(self.affine)
            self.dz = linepts[-1, 0] - linepts[0, 0] * (-1 / 10)
            self.dy = linepts[-1, 1] - linepts[0, 1] * (0.9434 / 10)
            self.dx = linepts[-1, 2] - linepts[0, 2] * (-1 / 10)
            self.vector = np.array([self.dz, self.dy, self.dx, 1])
            print('Affine Determinant', np.linalg.det(self.affine))
            det = np.linalg.det(self.affine)
            self.vectorPrime = self.affine.dot(self.vector)
            self.affScale = np.linalg.norm(self.vectorPrime[0:3]) / np.linalg.norm(self.vector)
            #self.affScale = 1
            print('Aff Scale', self.affScale)

            self.scale = 3840 / (self.trigScale * (10 / self.affScale))
            print('Scale', self.scale)
            self.scale /= 10
            channels_scale = np.linspace(0, self.scale / det, 384)
        
            newPoints = []
            for i in range(len(self.channelsOriginal)):
                newPoints.append([self.channelsOriginal[i][0], channels_scale[i]])
        
            self.ogPoints = newPoints
            self.channels = newPoints
            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
        else:
            newPoints = []
            for i in range(len(self.channelsOriginal)):
                newPoints.append([self.channelsOriginal[i][0], old_y[i][1]])

            self.ogPoints = newPoints
            self.channels = newPoints
            self.adj = [[i, i + 1] for i in range(len(self.channelsOriginal) - 1)]
            self.channelsPlot.setData(pos=np.array(self.channels, dtype=float), adj=np.array(self.adj, dtype=int))
            self.linearSpacePoints(points_added)

        """
        self.volArray = self.getColorVolume(intensity_values)
        print(self.volArray.shape)
        
        self.volArray *= 2
        rot = np.rot90(self.volArray)
        flip = np.flipud(rot)
      
        self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
        view = self.image.getView()

        self.points = [[100, t] for t in range(j)]
        self.plItem = pg.ScatterPlotItem(pos=self.points, pen=QtGui.QPen(QColor('red')), brush=QtGui.QBrush(QColor('red')))
        #self.plItem.setClickable(True)
        self.plItem.sigClicked.connect(self.onclickProbe)
        view.addItem(self.channelsPlot)
        view.addItem(self.plItem)
        view.addItem(self.textItem)
        """

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
        self.initUI()
    
    def initUI(self):
        self.setWindowTitle(self.title)
        self.width = int(4000) #default = coronal view
        self.height = int(4000)
        self.remove = True
        self.show = True
        self.prevProbe = ''

        im8 = np.ones((self.height,self.width),dtype='uint8')*255
        self.image = pg.image(im8)
        self.image.ui.histogram.hide()
        self.image.ui.roiBtn.hide()
        self.image.ui.menuBtn.hide()
        self.image.setObjectName('image')
        #self.image.setImage(im8.transpose())

        self.pointsAdded = []
        self.channelCoords = {}
        self.lineItems = []
        self.oldChannels = [] # stack for undoing lines

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
        self.volumeImage = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_red.mhd'))).T
        self.metricsList = ['Spread', 'Amplitude', 'Cumulative_Drift', 'Velocity_Above', 'Velocity_Below']

        self.unitPlot = PlotDisplayItem('unit_density', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID, self.metricsList)
        self.spreadPlot = PlotDisplayItem('spread', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID, self.metricsList)
        #self.firePlot = PlotDisplayItem('firing_rate', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID)
        self.velocityAbovePlot = PlotDisplayItem('velocity_above', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID, self.metricsList)
        self.velocityBelowPlot = PlotDisplayItem('velocity_below', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID, self.metricsList)
        self.ampPlot = PlotDisplayItem('amplitude', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID, self.metricsList)
        #self.repolarPlot = PlotDisplayItem('repolarization_slope', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID)
        #self.dPrimePlot = PlotDisplayItem('d_prime', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID)
        self.cumDriftPlot = PlotDisplayItem('cumulative_drift', self.waveform, self.metrics, self.volumeImage, self.probeAnnotations, self.mouseID, self.metricsList)
        self.plots = {'unit_density': self.unitPlot, 'spread': self.spreadPlot,'cumulative_drift': self.cumDriftPlot, 
                      'velocity_above': self.velocityAbovePlot, 'velocity_below': self.velocityBelowPlot, 'amplitude': self.ampPlot}
        
        for plot in self.plots:
            others = [k for k in self.plots if k != plot]
            for other in others:
                self.plots[plot].otherPlots.append(self.plots[other]) # add other plots to list for updating 

        #self.getColorVolume()

        # main layout with region of interest
        self.mainLayout = QVBoxLayout()
        self.imageLayout = QHBoxLayout()

        
        self.imageLayout.addWidget(self.image)

        self.probes = self.probeAnnotations['probe_name'].unique()
        self.probeLetters = [s[s.index(' ')+1:][0] for s in self.probes]
        self.probeLetters = sorted(list(set(self.probeLetters)))
        print(self.probeLetters)
        self.probeDropDown = QComboBox()
        for probe in sorted(self.probes):
            self.probeDropDown.addItem(probe)

        self.metrics = QComboBox()
        for metric in self.metricsList:
            self.metrics.addItem(metric)

        self.viewButton = QPushButton('View Probe Region with Selected Metric')
        self.viewButton.clicked.connect(self.displayRegion)

        self.toggleProbeButton = QPushButton('Toggle Probe')
        self.toggleProbeButton.clicked.connect(self.toggleProbe)

        self.resetPlotButton = QPushButton('Reset Metric Plot')
        self.resetPlotButton.clicked.connect(self.resetPlot)

        self.warpButton = QPushButton('Warp to CCF')
        self.warpButton.clicked.connect(self.warpChannels)

        self.probeViewLayout = QHBoxLayout()
        self.probeViewLayout.addWidget(self.probeDropDown)
        self.probeViewLayout.addWidget(self.metrics)
        self.probeViewLayout.addWidget(self.viewButton)
        #self.probeViewLayout.addWidget(self.toggleProbeButton)
        self.probeViewLayout.addWidget(self.resetPlotButton)

        self.probeViewLayout.addWidget(self.warpButton)
        self.probeViewLayout.setAlignment(QtCore.Qt.AlignTop)
        self.mainLayout.addLayout(self.imageLayout)
        self.mainLayout.addLayout(self.probeViewLayout)

        self.setLayout(self.mainLayout)
        self.showMaximized()


    def resetPlot(self):
        self.plots['unit_density'].resetPlot()
        self.plots[self.metrics.currentText().lower()].resetPlot()

        view = self.image.getView()
        for item in self.lineItems:
            view.removeItem(item)

        """
        if hasattr(self, 'plItem'):
            view.removeItem(self.plItem)
        """

        self.lineItems.clear()
        self.pointsAdded.clear()
        self.oldChannels.clear()

    def toggleProbe(self):
        print(self.show)
        view = self.image.getView()
        if self.show:
            view.removeItem(self.plItem)
            self.show = False
        else:
            view.addItem(self.plItem)
            self.show = True

    def getColorVolume(self, intensity_values, rgb_levels=DEFAULT_COLOR_VALUES):
       colarray = np.clip(intensity_values, a_min=rgb_levels[0][0], a_max=rgb_levels[0][1]) - rgb_levels[0][0]
       colarray = (colarray * 255. / (rgb_levels[0][1] - rgb_levels[0][0])).astype('uint8')
       
       return colarray

    # warps the channels to the ccf
    def warpChannels(self):
        self.channelsOriginal = self.unitPlot.channelsOriginal
        self.channels = self.unitPlot.channels
        #self.coords = self.unitPlot.coords

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
            channel_dict['DV'].append(coord[1] * (3840 / (100 * self.plots['unit_density'].affScale)))
            channel_dict['ML'].append(coord[2])

        probe_name = self.probeDropDown.currentText()
        probe = [probe_name for i in range(len(channel_dict['AP']))]
        channel_dict['probe_name'] = probe
        df_channel = pd.DataFrame(channel_dict)
        
        warp_channels(self.storageDirectory, df_channel, self.field, self.reference, self.probeDropDown.currentText(), self.mouseID, channels)

    # displays the region of interest surrounding the probe track
    def displayRegion(self):
        probe = self.probeDropDown.currentText()
        metric = self.metrics.currentText().lower()

        if self.prevProbe == '' or self.prevProbe != probe:
            self.updateDisplay(probe)
            probe_let_num = probe[probe.index(' ')+1:]
            days = sorted(list(self.waveMetricsPath.keys()))
            key = days[int(probe_let_num[1]) - 1]
            paths = self.waveMetricsPath[key]

            for p in paths:
                if 'probe' + probe_let_num[0] in p:
                    path = p


            #for plot in self.plots:
            self.plots['unit_density'].updateMetrics(path)
            #self.plots['unit_density'].resetPlot(remove_probe=True)
            self.plots['unit_density'].updateDisplay(probe, self.linepts, self.intensityValues)

            #for plot in self.plots:
            self.plots[metric].updateMetrics(path)
            #self.plots[metric].resetPlot(remove_probe=True)
            self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues)

            self.resetPlot()
            self.prevProbe = probe
            self.oldMetric = metric.lower()
        elif self.prevProbe == probe: # just update with new metric
            view = self.image.getView()
            view.removeItem(self.plots[self.oldMetric].channelsPlot)

            #self.updateDisplay(probe)
            probe_let_num = probe[probe.index(' ')+1:]
            days = sorted(list(self.waveMetricsPath.keys()))
            key = days[int(probe_let_num[1]) - 1]
            paths = self.waveMetricsPath[key]
            for p in paths:
                if 'probe' + probe_let_num[0] in p:
                    path = p

            self.plots[metric].updateMetrics(path)
            #self.plots[metric].resetPlot(remove_probe=True)
            self.plots[metric].updateDisplay(probe, self.linepts, self.intensityValues, keep_y=True, old_y=self.plots[self.oldMetric].channels, points_added=self.pointsAdded)
            self.oldMetric = metric.lower()
            self.updateDisplay(probe)

    def removeLineHelper(self, y_coord, ind):
        view = self.image.getView()
        item = self.lineItems[ind]
        view.removeItem(item)
        self.lineItems.remove(item)
        self.pointsAdded.remove(y_coord)

        """
        if len(self.pointsAdded) == 0:
            self.plots['unit_density'].channels = self.plots['unit_density'].ogPoints
            self.plots[self.metrics.currentText().lower()].channels = self.plots[self.metrics.currentText().lower()].ogPoints
        """
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

    # updates the other plots to align it with the plot that was clicked on
    def updateAfterClick(self, other_plot, line_point, channel, pointsAdded):
        print('Other plot clicked', other_plot.measurement)
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
                    else:
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

    def onClickProbeHelper(self, click_plot, line_point, is_unit=True):
        channel = click_plot.channels.index([click_plot.channelsPlot.scatterPoint[0], click_plot.channelsPlot.scatterPoint[1]])
        click_plot.oldChannels.append(click_plot.channels)
        flag = True

        if line_point[1] != click_plot.channelsPlot.scatterPoint[1]:
            if line_point[1] < click_plot.channelsPlot.scatterPoint[1]:
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
                    elif len(greater) > 0 and len(less) == 0: # new line is below exiting, so shift
                        for p in click_plot.channels:
                            if p[1] <= click_plot.channelsPlot.scatterPoint[1]:
                                newPoints.append([p[0], p[1] - diff])
                            else:
                                newPoints.append([p[0], p[1]])
                    else:
                        flag = click_plot.linspaceAbove(srt, channel, diff, flag)

                        for i in range(len(click_plot.channels)):
                            p = click_plot.channels[i]
                            if p[1] > click_plot.channelsPlot.scatterPoint[1]:
                                click_plot.channels[i] = [p[0], p[1] - diff]
            elif line_point[1] > click_plot.channelsPlot.scatterPoint[1]:
                diff = line_point[1] - click_plot.channelsPlot.scatterPoint[1]
                newPoints = []

                if len(self.pointsAdded) == 0:
                    for p in click_plot.channels:
                        newPoints.append([p[0], p[1] + diff])
                else:
                    srt = sorted(self.pointsAdded)
                    greater = [t for t in self.pointsAdded if t > line_point[1]]
                    less = [t for t in self.pointsAdded if t < line_point[1]]

                    if len(less) > 0 and len(greater) > 0:
                        flag = click_plot.linspaceBelow(srt, channel, diff, flag)
                    elif len(greater) > 0 and len(less) == 0:
                        flag = click_plot.linspaceBelow(srt, channel, diff, flag)

                        for i in range(len(click_plot.channels)):
                            p = click_plot.channels[i]
                            if p[1] < click_plot.channelsPlot.scatterPoint[1]:
                                click_plot.channels[i] = [p[0], p[1] + diff]
                    else:
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
        else:
            pts = [[t, line_point[1]] for t in range(int(click_plot.channelsPlot.scatterPoint[0]), int(line_point[0]))]
            diff = 0
            h_line = pg.ScatterPlotItem(pos=pts, pen=QtGui.QPen(QColor('yellow')), brush=QtGui.QBrush(QColor('yellow')), size=2)
            self.lineItems.append(h_line)

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

            for k in range(-80,80):
                try:
                    intensity_values[j,k+80] = (self.volumeImage[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                except IndexError:
                    pass
        
        # display image
        #intensity_values *= 2.
        #self.trigScale = intensity_values.shape[0] / intensity_values.shape[1] 
        self.volArray = self.getColorVolume(intensity_values)
        print(self.volArray.shape)
        
        self.volArray *= 2
        rot = np.rot90(self.volArray)
        flip = np.flipud(rot)
      
        self.image.setImage(flip[:, 100:], levels=(0, 255), autoRange=False)
        view = self.image.getView()

        self.points = [[100, t] for t in range(j)]
        self.plItem = pg.ScatterPlotItem(pos=self.points, pen=QtGui.QPen(QColor('red')), brush=QtGui.QBrush(QColor('red')))
        #self.plItem.setClickable(True)
        self.plItem.sigClicked.connect(self.onclickProbe)

        view.addItem(self.plots['unit_density'].channelsPlot)
        view.addItem(self.plots[self.metrics.currentText().lower()].channelsPlot)

        view.addItem(self.plItem)
        #view.addItem(self.textItem)

if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID

    app = QApplication(sys.argv)
    w = VolumeAlignment(mouse_id)
    sys.exit(app.exec_())