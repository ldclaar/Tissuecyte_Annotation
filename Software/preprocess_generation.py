import pandas as pd
import numpy as np
import argparse
import pathlib
import os
import SimpleITK as sitk
from sklearn.cluster import KMeans
from psycopg2 import connect, extras
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=False)
DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]

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

class generateImages():
    def __init__(self, mouse_id):
        self.mouseID = mouse_id
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_path)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')

        self.reference = sitk.ReadImage(self.reference_file)
        self.field = sitk.ReadImage(self.field_file)

        if os.path.exists(os.path.join(self.storageDirectory, 'reassigned', 'probe_annotations_{}_reassigned.csv'.format(self.mouseID))):
            self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'reassigned', 'probe_annotations_{}_reassigned.csv'.format(self.mouseID)))
        else:
            self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))

        with open(os.path.join(self.workingDirectory, 'field_reference', 'name_map.pkl'), 'rb') as f:
            self.name_map = pickle.load(f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'acrnm_map.pkl'), 'rb') as f:
            self.acrnm_map = pickle.load(f)
        with open(os.path.join(self.workingDirectory, 'field_reference', 'color_map.pkl'), 'rb') as f:
            self.colormap = pickle.load(f)
        self.anno = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.workingDirectory, 'field_reference', 'ccf_ano.mhd')))

        self.volumeRed = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_red.mhd'))).T
        self.volumeGreen = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_green.mhd'))).T
        self.volumeBlue = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_blue.mhd'))).T
        self.myFont = ImageFont.load_default()

    # helper function to normalize slice
    def getColorVolume(self, rgb_levels=DEFAULT_COLOR_VALUES):
        level_adjusted_arrays = []
        for colori, int_level in zip(['red', 'green', 'blue'], rgb_levels):
            colarray = np.clip(self.int_arrays[colori], a_min=int_level[0], a_max=int_level[1]) - int_level[0]
            colarray = (colarray * 255. / (int_level[1] - int_level[0])).astype('uint8')
            level_adjusted_arrays.append(colarray)
        return np.stack(level_adjusted_arrays, axis=-1)

    # generates images before alignment
    # probe: string, the probe to be displayed
    def imageGenerate(self, probe):
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
            
        print(linepts.shape[0])
        # extract region of image, build row cells containing image
        intensity_values_red = np.zeros((linepts.shape[0],160))
        intensity_values_green = np.zeros((linepts.shape[0],160))
        intensity_values_blue = np.zeros((linepts.shape[0],160))

        # for mask displaying regions
        intensity_ccf_red = np.zeros((linepts.shape[0],160))
        intensity_ccf_green = np.zeros((linepts.shape[0],160))
        intensity_ccf_blue = np.zeros((linepts.shape[0],160))
        #volume_warp = np.zeros((528, 320, 456))

        self.coords = {}
        self.linepts = linepts
        self.intensityValues = intensity_values_red
        coords_int = {}
        labels_pos = {}
        keys = list(self.acrnm_map.keys())
        values = list(self.acrnm_map.values())

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2])

            for k in range(-80,80): # build slice that probe cuts through
                try:
                    point = [int(np.round(linepts[j, 0] / 2.5)), int(np.round(linepts[j, 1] / 2.5)), int(np.round((linepts[j, 2] + k) / 2.5))]
                    
                    intensity_values_red[j,k+80] = (self.volumeRed[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    intensity_values_green[j,k+80] = (self.volumeGreen[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    intensity_values_blue[j,k+80] = (self.volumeBlue[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])

                    coords_int[(j, k+80)] = point
                    struct = self.anno[point[0], point[1], point[2]]

                    if struct in values and j >= 400:
                        ind = values.index(struct)
                        key = keys[ind]

                        if not key[0].islower():
                            if key in labels_pos:
                                labels_pos[key].append((j, k+80))
                            else:
                                labels_pos[key] = [(j, k+80)]
                except IndexError:
                    pass
        
        for row in range(linepts.shape[0]):
            for col in range(160):
                try:
                    if(row, col) in coords_int:
                        point = coords_int[(row, col)]
                        struct = self.anno[point[0], point[1], point[2]]

                        if struct in self.colormap:
                            intensity_ccf_red[row,col] = self.colormap[struct][0]
                            intensity_ccf_green[row,col] = self.colormap[struct][1]
                            intensity_ccf_blue[row,col] = self.colormap[struct][2]
                except IndexError:
                    pass

        # volume that probe cuts through
        self.int_arrays = {}
        self.int_arrays['red'] = intensity_values_red
        self.int_arrays['green'] = intensity_values_green
        self.int_arrays['blue'] = intensity_values_blue
        self.volArray = self.getColorVolume()
        print(self.volArray.shape)

        level_arrays = []
        level_arrays.append(intensity_ccf_red)
        level_arrays.append(intensity_ccf_green)
        level_arrays.append(intensity_ccf_blue)

        # mask of ccf regions 
        mask = np.stack(level_arrays, axis=-1)
        print(mask.shape)

        self.mask = Image.fromarray(mask.astype(np.uint8)).copy()
        im = Image.fromarray(self.volArray).copy()

        overlay = Image.blend(im, self.mask, 0.5).copy()
        draw = ImageDraw.Draw(overlay)
        draw_mask = ImageDraw.Draw(self.mask)

        for key in labels_pos:
            pos = np.array(labels_pos[key])
            center = pos.mean(axis=0)

            draw_mask.text((center[1], center[0]), key, fill=(255, 255, 255), font=self.myFont)
            draw.text((center[1], center[0]), key, fill=(255, 255, 255), font=self.myFont)

        im.save(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_slice.png'.format(probe.replace(' ', '_')))) # save slice
        self.mask.save(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_mask.png'.format(probe.replace(' ', '_')))) # save mask with text
        overlay.save(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_overlay.png'.format(probe.replace(' ', '_')))) # save overlay
        with open(os.path.join(self.workingDirectory, self.mouseID, 'images', '{}_labels.pickle'.format(probe.replace(' ', '_'))), 'wb') as handle:
            pickle.dump(labels_pos, handle, protocol=pickle.HIGHEST_PROTOCOL)
        """
        reference_dict = {'AP': [], 'DV': [], 'ML': []}
        mask_points = {}
        # warp 3d point at each pixel to ccf using deformation field
        for row in range(linepts.shape[0]):
            for col in range(intensity_values_red.shape[1]):
                try:
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
        """

if __name__ == '__main__':
    args = parser.parse_args()
    mouse_id = args.mouseID

    #mouse_ids = ['604914', '608672', '612090', '614547', '607660', '614608', '615047', '615048', '615564', '623322', '623784', '623319',
                 #'623786', '626279', '632296']
    mouse_ids = ['623786', '626279', '632296']

    for mid in mouse_ids:
        print('MouseID', mid)
        pre = generateImages(mid)
        probes = sorted(pre.probeAnnotations['probe_name'].unique())

        for probe in probes:
            pre.imageGenerate(probe)
        print()