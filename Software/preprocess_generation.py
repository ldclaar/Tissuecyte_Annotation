import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pathlib
import os
import SimpleITK as sitk
from sklearn.cluster import KMeans
from psycopg2 import connect, extras

# query lims and return result from given query
def query_lims(query_string):
    con = connect(
    dbname='lims2',
    user='limsreader',
    host='limsdb2',
    password='limsro',
    port=5432,
    )
    con.set_session(
        readonly=True, 
        autocommit=True,
    )
    cursor = con.cursor(
        cursor_factory=extras.RealDictCursor,
    )
    cursor.execute(query_string)
    result = cursor.fetchall()

    return result

# gets the tissuecyte info for the mouse id
def get_tc_info(mouse_id):
    TISSUECYTE_QRY = '''
            SELECT *
            FROM image_series im
            WHERE im.specimen_id = {}
        '''

    storage_directory = ''
    tc = query_lims(TISSUECYTE_QRY.format(get_specimen_id_from_labtracks_id(int(mouse_id))))
    for row in tc:
        d = dict(row)
        if d['alignment3d_id'] != None:
            storage_directory = d['storage_directory']

    print(storage_directory)
    return storage_directory

def get_specimen_id_from_labtracks_id(labtracks_id):
    SPECIMEN_QRY = '''
            SELECT *
            FROM specimens sp
            WHERE sp.external_specimen_name=cast({} as character varying)
        '''

    mouse_info = query_lims(SPECIMEN_QRY.format(int(labtracks_id)))
    return mouse_info[0]['id']

# cluster based on number of annotations per probe
def cluster_annotations(k, nonzero_points, warped_dict, all_dict):
    #print('K', k)
    kmeans = KMeans(n_clusters=k, random_state=0)
    labels = kmeans.fit_predict(nonzero_points)
    #print(labels)
    centers = kmeans.cluster_centers_

    for center in centers:
        ap = int(np.round(center[2]))
        dv = int(np.round(center[1]))
        ml = int(np.round(center[0]))

        warped_dict['AP'].append(ap)
        warped_dict['DV'].append(dv)
        warped_dict['ML'].append(ml)

        if all_dict != None:
            all_dict['AP'].append(ap)
            all_dict['DV'].append(dv)
            all_dict['ML'].append(ml)

# warps the volume using the reference, field and interpolator
def warp_execute( volume, reference, field, interpolator ) :
    deformWarped = sitk.Warp( volume, field, sitk.sitkLinear, 
                        reference.GetSize(), reference.GetOrigin(), reference.GetSpacing())
    return deformWarped

class generateImages():
    def __init__(self, mouse_id):
        self.mouseID = mouse_id
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_file)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')

        self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))

    # displays the region along the probe track
    # probe: string, the probe to be displayed from the drop down
    def updateDisplay(self, probe):
        self.showMask = True
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

        self.coords = {}
        self.linepts = linepts
        coords_int = {}

        print(linepts.shape[0])
        # extract region of image, build row cells containing image
        intensity_values_red = np.zeros((linepts.shape[0],160))
        intensity_values_green = np.zeros((linepts.shape[0],160))
        intensity_values_blue = np.zeros((linepts.shape[0],160))
        self.intensityValues = intensity_values_red

        # for mask displaying regions
        intensity_ccf_red = np.zeros((linepts.shape[0],160))
        intensity_ccf_green = np.zeros((linepts.shape[0],160))
        intensity_ccf_blue = np.zeros((linepts.shape[0],160))
        volume_warp = np.zeros((528, 320, 456))

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2])

            for k in range(-80,80): # build slice that probe cuts through
                try:
                    ccf_point = [int(np.round(linepts[j, 0] / 2.5)), int(np.round(linepts[j, 1] / 2.5)), int(np.round((linepts[j, 2] + k) / 2.5))]
                    intensity_values_red[j,k+80] = (self.volumeRed[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    intensity_values_green[j,k+80] = (self.volumeGreen[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    intensity_values_blue[j,k+80] = (self.volumeBlue[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])

                    coords_int[(j, k+80)] = ccf_point
                except IndexError:
                    pass

        reference_dict = {'AP': [], 'DV': [], 'ML': []}
        mask_points = {}
        # warp 3d point at each pixel to ccf using deformation field
        try:
            for row in range(linepts.shape[0]):
                for col in range(intensity_values_red.shape[1]):
                    point = coords_int[row, col]
                    volume_warp[point[0], point[1], point[2]] = 1

                    vol_im = sitk.GetImageFromArray(volume_warp.T)
                    vol_im.SetSpacing((25, 25, 25))

                    warped_im = warp_execute(vol_im, self.reference, self.field, sitk.sitkLinear)
                    warped_arr = sitk.GetArrayFromImage(warped_im)
                    nonzero_points = np.argwhere(warped_arr > 0)
                    cluster_annotations(1, nonzero_points, reference_dict, None)
                    mask_points[(row, col)] = [reference_dict['AP'][-1], reference_dict['DV'][-1], reference_dict['ML'][-1]]
                    volume_warp[volume_warp > 0] = 0 # reset array
        except IndexError:
            pass

if __name__ == '__main__':
    pre = generateImages()
    probes = pre.probeAnnotations['probe_name'].unique()

    pre.updateDisplay(probes[0])