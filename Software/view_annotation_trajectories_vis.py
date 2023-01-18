import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, QAbstractTableModel
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, QTableView, \
    QVBoxLayout, QWidget, QPushButton, QGridLayout, QCheckBox, QFormLayout, QFileDialog, QComboBox, QMessageBox, QLineEdit
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtGui import QColor
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import sys
import pyqtgraph as pg
import os
from matplotlib import cm
import argparse
import pathlib
import SimpleITK as sitk
from sklearn.cluster import KMeans
import visvis as vis
#class pandasModel(QAbstractTableModel):

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)

backend = 'pyqt5'
app = vis.use(backend)

class AnnotationProbesViewer(QWidget):
    # initialize fields
    def __init__(self, mouse_id):
        super().__init__()
        # directory and csv fields
        self.mouseID = mouse_id
        
        self.dir = '//allen/programs/mindscope/workgroups/np-behavior/tissuecyte'
        self.workingDirectory = pathlib.Path('{}/{}'.format(self.dir, self.mouseID))
        print('Fetching Data')
        self.storage_directory = pathlib.Path('{}/field_reference'.format(self.dir))
        #model_directory = '//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56'
        self.model_directory = pathlib.Path('{}/field_reference'.format(self.dir))
        #input_file = os.path.join( working_directory, 'output_8', 'affine_vol_10.mhd')
    
        #affineAligned = sitk.ReadImage( input_file )
        #affineAligned.SetSpacing((25, 25, 25))
        #self.field_file = os.path.join( self.storage_directory, 'dfmfld.mhd')
        self.reference_file = os.path.join( self.model_directory, 'average_template_25.nrrd')

        self.volume = sitk.GetArrayFromImage(sitk.ReadImage(self.reference_file))
        self.reference = sitk.ReadImage( self.reference_file )
        #self.field = sitk.ReadImage( self.field_file )

        self.annotations = pd.read_csv(os.path.join(self.workingDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))
        self.updatedAnnotations = self.annotations.copy(deep=True)

        self.rgb = {'A1': '(255, 228, 225)', 'A2': '(255, 0, 0)', 'A3': '(240, 128, 128)', 'A4': '(139, 0, 0)',
                    'B1': '(173, 216, 230)', 'B2': '(0, 0, 255)', 'B3': '(70, 130, 180)', 'B4': '(0, 0, 139)',
                    'C1': '(255, 192, 203)', 'C2': '(255, 0, 255)', 'C3': '(218, 112, 214)', 'C4': '(255, 20, 147)',
                    'D1': '(210, 180, 140)', 'D2': '(255, 128, 0)', 'D3': '(255, 215, 0)', 'D4': '(218, 165, 32)',
                    'E1': '( 0, 255, 255)', 'E2': '(95, 158, 160)', 'E3': '(127, 255, 212)', 'E4': '(72, 209, 204)',
                    'F1': '(144, 238, 144)', 'F2': '(0, 128, 0)', 'F3': '(128, 128, 0)', 'F4': '(107, 142, 35)'}

        self.probe_lines = {}
        # Make figure using "self" as a parent
        Figure = app.GetFigureClass()
        self.fig = Figure(self)
        print(type(self.fig._widget))
        #self.fig.bgcolor = 'k'
        #self.panel = QWidget(self)
        self.create_main_frame()
        self.ax = vis.subplot(111)
        # layouts for the display
        self.labelLayout = QVBoxLayout()
        self.labelLayout.addWidget(self.main_frame)
        self.labelOld = QHBoxLayout()
        self.labelNew = QHBoxLayout()
        self.scatterLayout = QVBoxLayout()
        self.mainLayout = QHBoxLayout()
        self.hButtonLayout = QHBoxLayout()
        #self.fig._widget.setAutoFillBackground(True)
        
        #self.mainLayout.addidget(self.panel, 1)
        self.mainLayout.addWidget(self.fig._widget)

        # display for current probe and number drop down
        self.probeOld = QComboBox()
        self.probeOld.addItem('Current Probe')
        self.probeOld.addItem('A')
        self.probeOld.addItem('B')
        self.probeOld.addItem('C')
        self.probeOld.addItem('D')
        self.probeOld.addItem('E')
        self.probeOld.addItem('F')
        self.labelOld.addWidget(self.probeOld)

        self.trialOld = QComboBox()
        self.trialOld.addItem('Current Number')
        self.trialOld.addItem('1')
        self.trialOld.addItem('2')
        self.trialOld.addItem('3')
        self.trialOld.addItem('4')
        self.trialOld.addItem('5')
        self.labelOld.addWidget(self.trialOld)
        
        # display for new probe and number drop down
        self.probeNew = QComboBox()
        self.probeNew.addItem('New Probe')
        self.probeNew.addItem('A')
        self.probeNew.addItem('B')
        self.probeNew.addItem('C')
        self.probeNew.addItem('D')
        self.probeNew.addItem('E')
        self.probeNew.addItem('F')
        self.labelNew.addWidget(self.probeNew)

        self.trialNew = QComboBox()
        self.trialNew.addItem('New Number')
        self.trialNew.addItem('1')
        self.trialNew.addItem('2')
        self.trialNew.addItem('3')
        self.trialNew.addItem('4')
        self.labelNew.addWidget(self.trialNew)

        self.switchButton = QPushButton('Switch Probes')
        self.switchButton.clicked.connect(self.switchProbes)


        self.updateButton = QPushButton('Reassign Probe')
        self.updateButton.clicked.connect(self.updateProbe)
        self.labelLayout.addLayout(self.labelOld)
        self.labelLayout.addLayout(self.labelNew)

        self.reassignSwitchLayout = QHBoxLayout()
        self.reassignSwitchLayout.addWidget(self.updateButton)
        self.reassignSwitchLayout.addWidget(self.switchButton)
        #self.labelLayout.addWidget(self.updateButton)

        self.checkLabel = QLabel('Select probes to show/hide')
        self.checkLayout = QHBoxLayout()
        
        self.aCheckBox = QCheckBox('A Probes')
        self.bCheckBox = QCheckBox('B Probes')
        self.cCheckBox = QCheckBox('C Probes')
        self.dCheckBox = QCheckBox('D Probes')
        self.eCheckBox = QCheckBox('E Probes')
        self.fCheckBox = QCheckBox('F Probes')

        self.checkLayout.addWidget(self.aCheckBox)
        self.checkLayout.addWidget(self.bCheckBox)
        self.checkLayout.addWidget(self.cCheckBox)
        self.checkLayout.addWidget(self.dCheckBox)
        self.checkLayout.addWidget(self.eCheckBox)
        self.checkLayout.addWidget(self.fCheckBox)

        self.probe_checks = {'A': self.aCheckBox, 'B': self.bCheckBox, 'C': self.cCheckBox, 'D': self.dCheckBox, 'E': self.eCheckBox, 'F': self.fCheckBox}
        self.toggleButton = QPushButton('Update Probe Trajectory Display')
        self.toggleButton.clicked.connect(self.updateDisplay)

        self.saveButton = QPushButton('Save Updated Probes')
        self.saveButton.clicked.connect(self.saveUpdatedProbes)


        """
        # toggle layouts for each probe
        self.aCheckLayout = QGridLayout()
        self.bCheckLayout = QGridLayout()
        self.cCheckLayout = QGridLayout()
        self.dCheckLayout = QGridLayout()
        self.eCheckLayout = QGridLayout()
        self.fCheckLayout = QGridLayout()
        self.checkLayout = QVBoxLayout()
        self.probe_checks = {} # dictionary to see if probe has been checked or not

        for i in range(1, 7):
            a_cbox = QCheckBox('A%d' %(i))
            b_cbox = QCheckBox('B%d' %(i))
            c_cbox = QCheckBox('C%d' %(i))
            d_cbox = QCheckBox('D%d' %(i))
            e_cbox = QCheckBox('E%d' %(i))
            f_cbox = QCheckBox('F%d' %(i))

            self.aCheckLayout.addWidget(a_cbox, 0, i)
            self.bCheckLayout.addWidget(b_cbox, 0, i)
            self.cCheckLayout.addWidget(c_cbox, 0, i)
            self.dCheckLayout.addWidget(d_cbox, 0, i)
            self.eCheckLayout.addWidget(e_cbox, 0, i)
            self.fCheckLayout.addWidget(f_cbox, 0, i)

            self.probe_checks['Probe A%d' %(i)] = a_cbox
            self.probe_checks['Probe B%d' %(i)] = b_cbox
            self.probe_checks['Probe C%d' %(i)] = c_cbox
            self.probe_checks['Probe D%d' %(i)] = d_cbox
            self.probe_checks['Probe E%d' %(i)] = e_cbox
            self.probe_checks['Probe F%d' %(i)] = f_cbox

        self.checkLabel = QLabel('Select probes to show/hide')

        self.checkLayout.addWidget(self.checkLabel)
        self.checkLayout.addLayout(self.aCheckLayout)
        self.checkLayout.addLayout(self.bCheckLayout)
        self.checkLayout.addLayout(self.cCheckLayout)
        self.checkLayout.addLayout(self.dCheckLayout)
        self.checkLayout.addLayout(self.eCheckLayout)
        self.checkLayout.addLayout(self.fCheckLayout)
        self.toggleButton = QPushButton('Update Probe Trajectory Display')
        self.toggleButton.clicked.connect(self.updateDisplay)
        self.checkLayout.addWidget(self.toggleButton)
        """
        self.labelLayout.addLayout(self.reassignSwitchLayout)
        self.labelLayout.addWidget(self.checkLabel)
        self.labelLayout.addLayout(self.checkLayout)
        #self.labelLayout.addWidget(self.toggleButton)
        self.hButtonLayout.addWidget(self.toggleButton)
        self.hButtonLayout.addWidget(self.saveButton)
        self.labelLayout.addLayout(self.hButtonLayout)
        self.labelLayout.setAlignment(QtCore.Qt.AlignBottom)

        self.scatterLayout.addWidget(self.main_frame)
        self.scatterLayout.addLayout(self.labelLayout)
        
        self.qProbe = []
        self.qTrial = []

        self.update_plot(self.annotations)
        self.update_plot_2d(self.annotations)
        #self.mainLayout.addWidget(self.main_frame)
        self.mainLayout.addLayout(self.scatterLayout)
        self.setLayout(self.mainLayout)
        self.showMaximized()

    def switchProbes(self):
        old_probe = 'Probe' + ' ' + self.probeOld.currentText() + self.trialOld.currentText()
        new_probe = 'Probe' + ' ' + self.probeNew.currentText() + self.trialNew.currentText()

        self.probeOld.setCurrentText('Current Probe')
        self.trialOld.setCurrentText('Current Number')
        self.probeNew.setCurrentText('New Probe')
        self.trialNew.setCurrentText('New Number')

        self.updatedAnnotations.loc[self.updatedAnnotations['probe_name'] == new_probe, 'probe_name'] = 'temp'
        self.updatedAnnotations.loc[self.updatedAnnotations['probe_name'] == old_probe, 'probe_name'] = new_probe
        self.updatedAnnotations.loc[self.updatedAnnotations['probe_name'] == 'temp', 'probe_name'] = old_probe

        self.update_plot(self.updatedAnnotations)
        self.update_plot_2d(self.updatedAnnotations)

    def saveUpdatedProbes(self):
        if not os.path.exists(os.path.join(self.workingDirectory, 'reassigned')):
            os.mkdir(os.path.join(self.workingDirectory, 'reassigned'))

        self.updatedAnnotations.to_csv(os.path.join(self.workingDirectory, 'reassigned', 'probe_annotations_{}_reassigned.csv'.format(self.mouseID)))

        probes = self.updatedAnnotations['probe_name'].unqiue()

        for probe in probes:
            df_probe = self.updatedAnnotations.loc[self.updatedAnnotations['probe_name'] == probe]
            df_probe.to_csv(os.path.join(self.workingDirectory, 'reassigned', '{}_annotations_{}_reassigned.csv'.format(probe, self.mouseID)))

    # populates the drop down, used when probes are updated
    def populateDropDown(self, qProbe, qTrial):
        self.probeOld.clear()
        self.trialOld.clear()

        print('QProbe', qProbe)
        print('QTrial', qTrial)

        self.probeOld.addItem('Current Probe')
        self.probeOld.addItem('A')
        self.probeOld.addItem('B')
        self.probeOld.addItem('C')
        self.probeOld.addItem('D')
        self.probeOld.addItem('E')
        self.probeOld.addItem('F')

        self.trialOld.addItem('Current Number')
        self.trialOld.addItem('1')
        self.trialOld.addItem('2')
        self.trialOld.addItem('3')
        self.trialOld.addItem('4')

        for probe in qProbe: 
            self.probeOld.addItem(probe)
          
        for trial in qTrial:
            self.trialOld.addItem(trial)

    def updateDisplay(self):
        self.update_plot(self.updatedAnnotations)
        self.update_plot_2d(self.updatedAnnotations)

    # called when probes need to be updated
    def updateProbe(self):
        probes = self.updatedAnnotations['probe_name'].unique()
        print(probes)
        orig_counts = len([p for p in probes if 'orig_Probe' in p])
        
        if not 'orig_Probe' in self.probeOld.currentText(): # new probe is not one of the probes currently displayed
            probe_name_old = 'Probe' + ' ' + self.probeOld.currentText() + self.trialOld.currentText()
            probe_name_new = 'Probe' + ' ' + self.probeNew.currentText() + self.trialNew.currentText()
        else:
            probe_name_old = self.probeOld.currentText() + ' ' + self.trialOld.currentText()
            self.qProbe.pop(0)
            self.qTrial.pop(0)
            print(self.qProbe)
            probe_name_new = 'Probe' + ' ' + self.probeNew.currentText() + self.trialNew.currentText()


        if probe_name_new in probes and probe_name_new[probe_name_new.index(' ') + 1:] != self.trialOld.currentText(): # need to add this to drop down, so it can be reassigned later on
            print(self.trialOld.currentText())
            #self.probeOld.addItem('orig_Probe')
            self.qProbe.append('orig_Probe')
            #self.trialOld.addItem(probe_name_new[probe_name_new.index(' ') + 1:])
            self.qTrial.append(probe_name_new[probe_name_new.index(' ') + 1:])
            self.updatedAnnotations.loc[self.updatedAnnotations['probe_name'] == probe_name_new, 'probe_name'] = 'orig_{}'.format(probe_name_new)
       
        self.populateDropDown(self.qProbe, self.qTrial)
        print(probe_name_old)
        self.updatedAnnotations.loc[self.updatedAnnotations['probe_name'] == probe_name_old, 'probe_name'] = probe_name_new
        self.update_plot(self.updatedAnnotations)
        self.update_plot_2d(self.updatedAnnotations)

        self.probeOld.setCurrentText('Current Probe')
        self.trialOld.setCurrentText('Current Number')
        self.probeNew.setCurrentText('New Probe')
        self.trialNew.setCurrentText('New Number')
        """
        print(probe_name_old, probe_name_new)
        self.updatedAnnotations = self.annotations.replace(probe_name_old, probe_name_new)
        self.updatedAnnotations = self.updatedAnnotations[['AP', 'ML', 'DV', 'probe_name']]
        self.updatedAnnotations.to_csv(os.path.join(self.directory, 'post_probe_annotations_updated.csv'))
        self.update_plot(self.updatedAnnotations)
        """
    
    # removes point when point is clicked on
    def removePoint(self, event):
        x_ml = int(np.round(event.xdata * 2.5)) # x coordinate of click
        y_dv = int(np.round(event.ydata * 2.5)) # y coordinate of click

        ml_data = self.updatedAnnotations['ML']
        dv_data = self.updatedAnnotations['DV']

        x_data = (np.abs(ml_data - x_ml))
        y_data = (np.abs(dv_data - y_dv))

        min_x_index = -1

        x_sorted = sorted(list(x_data.nsmallest(10).index))
        y_sorted = sorted(list(y_data.nsmallest(10).index))

        for ind in x_sorted:
            if ind in y_sorted:
                min_x_index = ind
                break

        self.updatedAnnotations.drop([min_x_index], inplace=True)
        self.updatedAnnotations.reset_index(drop=True, inplace=True)
        self.update_plot(self.updatedAnnotations) 
        self.update_plot_2d(self.updatedAnnotations)
        """
        if x_data[min_x_index] < y_data[min_y_index]:
            self.updatedAnnotations.drop([min_x_index], inplace=True)
            self.updatedAnnotations.reset_index(drop=True, inplace=True)
            self.update_plot(self.updatedAnnotations) 
            self.update_plot_2d(self.updatedAnnotations)
        else:
            self.updatedAnnotations.drop([min_y_index], inplace=True)
            self.updatedAnnotations.reset_index(drop=True, inplace=True)
            self.update_plot(self.updatedAnnotations) 
            self.update_plot_2d(self.updatedAnnotations)
        """
      
    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        
        self.dpi = 100
        self.fig_2d = Figure((10, 5), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig_2d)
        self.canvas.setParent(self.main_frame)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig_2d.add_subplot(111)
        # Bind the 'pick' event for clicking on one of the bars
        #
        self.fig_2d.canvas.mpl_connect('button_press_event', self.removePoint)
        #self.fig.canvas.mpl_connect('key_press_event', self.moveSlice)

        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

     # updates the display
    def update_plot_2d(self, probe_annotations):
        self.axes.clear()

        self.colors = {'A1': 'mistyrose', 'A2': 'red', 'A3': 'light coral', 'A4': 'dark red',
                       'B1': 'light blue', 'B2': 'blue', 'B3': 'steel blue', 'B4': 'dark blue',
                       'C1': 'pink', 'C2': 'magenta', 'C3': 'orchid', 'C4': 'deep pink',
                       'D1': 'tan', 'D2': 'orange', 'D3': 'gold', 'D4': 'goldenrod',
                       'E1': 'cyan', 'E2': 'cadet blue', 'E3': 'aquamarine', 'E4': 'medium turquoise',
                       'F1': 'light green', 'F2': 'green', 'F3': 'olive', 'F4': 'olive drab'}


        probes = probe_annotations['probe_name'].unique()
        #self.axes.contourf(self.reference[:, 0], self.reference[:, 1], self.image[:, 2])
        #self.axes.set_xticklabels([])
        #self.axes.set_yticklabels([])
        #self.axes.set_zticklabels([])
        self.axes.set_facecolor('lightgrey') 
        
        print(self.volume.shape)
        
        for probe_idx, probe_trial in enumerate(probes):
            probe = probe_trial[probe_trial.index(' ') + 1:]
            x = probe_annotations[probe_annotations.probe_name == probe_trial].ML / 2.5
            y = probe_annotations[probe_annotations.probe_name == probe_trial].DV / 2.5
        
            z = probe_annotations[probe_annotations.probe_name == probe_trial].AP / 2.5

            # get trajectory
            if len(z) > 0:
                data = np.vstack((z,y,x)).T
                datamean = data.mean(axis=0)
                D = data - datamean
                m1 = np.min(D[:,1]) * 2
                m2 = np.max(D[:,1]) * 2
                uu,dd,vv = np.linalg.svd(D)

                linepts = vv[0] * np.mgrid[-200:200:0.7][:,np.newaxis]
                linepts += datamean
            
                if linepts[-1,1] - linepts[0,1] < 0:
                    linepts = np.flipud(linepts)

            if probe[0] in self.probe_checks and 'orig_Probe' not in probe_trial:
                if not self.probe_checks[probe[0]].isChecked(): # display probe
                    self.axes.scatter(x,y,c=self.colors[probe].replace(' ', ''), s=15, alpha=0.95)
                    self.axes.plot(linepts[:,2],linepts[:,1],color=self.colors[probe].replace(' ', ''), alpha=0.5, label=probe_trial)
            else: # make probe light grey, needs to be reassigned later
                self.axes.scatter(x,y,c='lightgrey', s=5, alpha=0.95)
                self.axes.plot(linepts[:,2],linepts[:,1],color='lightgrey', alpha=0.5, label=probe_trial)
        
        #self.axes.legend(loc='upper left')
        #self.axes.set_xlabel('ML')
        #self.axes.set_ylabel('DV')
        self.axes.set_title('2D View of Probe Annotations')
        #self.axes.set_zlabel('AP')
        self.axes.imshow(self.volume[:, :, :].T.sum(axis=0), cmap='gray')
        #self.axes.plot_trisurf(self.volume[:319, 0], self.volume[:, 1], self.volume[:, 2])
        self.canvas.draw()
        plt.show()

    # updates the display
    def update_plot(self, probe_annotations):
        #vis.clf()
        #self.axes.clear()
        self.ax.bgcolor='k'
        self.ax.axis.axisColor = 'w'
        vis.cla()

        vis.volshow3(self.volume)
        probes = probe_annotations['probe_name'].unique()
        #self.axes.contourf(self.reference[:, 0], self.reference[:, 1], self.image[:, 2])
        #self.axes.set_xticklabels([])
        #self.axes.set_yticklabels([])
        #self.axes.set_zticklabels([])
        #self.axes.set_facecolor('black') 
        
        print(self.volume.shape)
        legend = []

        for probe_idx, probe_trial in enumerate(probes):
            probe = probe_trial[probe_trial.index(' ') + 1:]
            x = probe_annotations[probe_annotations.probe_name == probe_trial].ML / 2.5
            y = probe_annotations[probe_annotations.probe_name == probe_trial].DV / 2.5
        
            z = probe_annotations[probe_annotations.probe_name == probe_trial].AP / 2.5

            # get trajectory
            if len(z) > 0:
                data = np.vstack((z,y,x)).T
                datamean = data.mean(axis=0)
                D = data - datamean
                m1 = np.min(D[:,1]) * 2
                m2 = np.max(D[:,1]) * 2
                uu,dd,vv = np.linalg.svd(D)

                linepts = vv[0] * np.mgrid[-200:200:0.7][:,np.newaxis]
                linepts += datamean
            
                if linepts[-1,1] - linepts[0,1] < 0:
                    linepts = np.flipud(linepts)


                if probe[0] in self.probe_checks and 'orig_Probe' not in probe_trial:
                    if not self.probe_checks[probe[0]].isChecked(): # display probe
                        #self.axes.scatter(x,y,c=self.colors[probe].replace(' ', ''), s=5, alpha=0.95)
                        col = eval(self.rgb[probe])
                        color = tuple(t / 255 for t in col)
                        #vis.plot(x, y, z, mc=color, mw=5, ms='s', lw=0, mec=color, axes=self.ax)
                        #print(l.points)
                        vis.plot(linepts[:, 2], linepts[:, 1], linepts[:, 0], lw=3, lc=color, axes=self.ax)
                        legend.append(probe_trial)
                else: # make probe light grey, needs to be reassigned later
                    legend.append(probe_trial)
                    color_grey = (211, 211, 211)
                    color = tuple(t / 255 for t in color_grey)
                    vis.plot(linepts[:, 2], linepts[:, 1], linepts[:, 0], lw=3, lc=color_grey, axes=self.ax)
        
                self.probe_lines[probe_trial] = linepts

        self.ax.legend = legend
        #self.ax.legend.bgcolor = 'grey'
        self.ax.axis.xLabel = 'ML'
        self.ax.axis.yLabel = 'DV'
        self.ax.axis.zLabel = 'AP'

if __name__ == '__main__':
    args = parser.parse_args()
    mouse_id = args.mouseID

    app.Create()
    m = AnnotationProbesViewer(mouse_id)
    app.Run()
    