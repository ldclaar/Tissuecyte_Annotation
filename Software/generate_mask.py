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
import scipy.ndimage as ndi
from PIL import Image

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session')

DEFAULT_COLOR_VALUES = [[0, 3000], [0, 3000], [0, 1000]]
class volumeMask:
    def __init__(self, mouse_id):
        self.mouseID = mouse_id
        # important directories
        self.workingDirectory = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
        self.modelDirectory = pathlib.Path('//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56/atlasVolume')
        self.storageDirectory = pathlib.Path(os.path.join(self.workingDirectory, self.mouseID))

        self.mcc = MouseConnectivityCache(resolution=25)
        self.anno, self.meta = self.mcc.get_annotation_volume()
        self.rsp = self.mcc.get_reference_space()
        self.rsp.remove_unassigned(); # This removes ids that are not in this particular reference space
        self.name_map = self.rsp.structure_tree.get_name_map() # dictionary mapping ids to structure names
        self.acrnm_map = self.rsp.structure_tree.get_id_acronym_map() # dictionary mapping acronyms to ids
        self.colormap = self.rsp.structure_tree.get_colormap() # the colormap used for the allen 3D mouse atlas ontology

        self.field_path = os.path.join(get_tc_info(self.mouseID), 'local_alignment', 'dfmfld.mhd') # get the deformation field for the given mouse
        self.field_file = pathlib.Path('/{}'.format(self.field_path))
        print(self.field_file)
        self.reference_file = os.path.join( self.modelDirectory, 'average_template_25.nrrd')

        self.reference = sitk.ReadImage( self.reference_file )
        self.field = sitk.ReadImage( self.field_file, sitk.sitkVectorFloat64 )
        self.dfmfld_transform = sitk.DisplacementFieldTransform(self.field)

        self.probeAnnotations = pd.read_csv(os.path.join(self.storageDirectory, 'probe_annotations_{}.csv'.format(self.mouseID)))
        self.volumeRed = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_red.mhd'))).T
        self.volumeGreen = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_green.mhd'))).T
        self.volumeBlue = sitk.GetArrayFromImage(sitk.ReadImage(os.path.join(self.storageDirectory, 'resampled_blue.mhd'))).T
        self.regions_of_interest = self.get_structure_ids(['CTX'])

        self.generateRoot()
        self.updateDisplay('Probe B2')
        self.generateMask()

    def generateMask(self):
        area_image = self.volArray.copy()
        boundary_image = self.volArray.copy()
        twoD = 2.857
        size = 4.464

        for num,struct in enumerate(self.structures_in_section[self.sorted_indx]):
            if struct != 997: # There is no need to map the root structure here, since this was calculated previously
                # For each structure, the points (voxels) within the volume are acquired.
                mask = self.rsp.make_structure_mask([struct])
                points = np.array(np.where(mask)).T.astype(float)
        
                # These points are transformed using the deformation and affine transforms.
                dfp = np.array(list(map(self.dfmfld_transform.TransformPoint,points*25)))
                #tfp = np.array(list(map(aff_trans.TransformPoint,dfp)))
        
                # Since these points will map to several sections, only the points in our desired section are selected
                sections = np.round(dfp[:,2]/100)
                loc = np.where(sections == 10)
    
                # Map points to 2D section
                mapped_points = dfp[loc,:2]*twoD/(2**4)
                mapped_points = mapped_points.squeeze()
                next_points = mapped_points + size
                rounded_mapped_points,rounded_next_points = np.round(mapped_points).reshape((-1,2)),np.round(next_points).reshape((-1,2))
    
                # Place points on blank array and fill in gaps.
                blank = np.zeros((self.volArray.shape[0],self.volArray.shape[1]))
                for i,rounded_mapped_point in enumerate(rounded_mapped_points):
                    rounded_next_point = rounded_next_points[i]
                    rounded_mapped_point = rounded_mapped_point.astype(int)
                    rounded_next_point = rounded_next_point.astype(int)
                    blank[rounded_mapped_point[1]:rounded_next_point[1],rounded_mapped_point[0]:rounded_next_point[0]] = 1
                area = ndi.binary_closing(ndi.binary_fill_holes(blank).astype(int)).astype(int)
    
                # Place structre in area image
                color = self.colormap[struct]
                area_mask = np.where(area==1)
                for i in np.arange(3):
                    area_image[:,:,i][area_mask] = color[i]
    
                # Draw region boundaries in boundary image
                if struct in self.regions_of_interest:
                    inner = area.copy()
                    for i in np.arange(6):
                        inner = ndi.binary_erosion(inner)
                    boundary = area - inner.astype(int)
                    boundary_mask = np.where(boundary == 1)
                    for i in np.arange(3):
                        boundary_image[:,:,i][boundary_mask] = 0
                    boundary_image[:,:,2][boundary_mask] = 255

        # Save area and boundary images        
        area_img = Image.fromarray(area_image.astype(np.uint8))
        area_img.save('C:/Users/arjun.sridhar/area_{}.png'.format(self.mouseID))
        boundary_img = Image.fromarray(boundary_image.astype(np.uint8))
        boundary_img.save('C:/Users/arjun.sridhar/boundary_{}.png'.format(self.mouseID))

        # Create and save area composite
        cbg = Image.open('C:/Users/arjun.sridhar/boundary_{}.png'.format(self.mouseID))
        fg = Image.open('C:/Users/arjun.sridhar/area_{}.png'.format(self.mouseID))
        cbg_alpha = Image.new("L",cbg.size,255)
        fg_alpha = Image.new("L",cbg.size,80)
        cbg.putalpha(cbg_alpha)
        fg.putalpha(fg_alpha)
        comp = Image.alpha_composite(cbg,fg)
        comp.save('C:/Users/arjun.sridhar/color_compositie_{}.png'.format(self.mouseID))

    # generates the root structure from the allen sdk
    def generateRoot(self):
        root_mask = self.rsp.make_structure_mask([997])
        root_points = np.array(np.where(root_mask)).T.astype(float)
        dfp = np.array(list(map(self.dfmfld_transform.TransformPoint,root_points*25)))
        sections = np.round(dfp[:,2]/100)
        loc = np.where(sections == 62)
        # Find annotated ids in these slices
        anno_loc = np.unique(root_points[loc][:,0])
        self.structures_in_section = np.unique(self.anno[anno_loc.astype(int),:,:])[1:]
        
        for region in self.regions_of_interest:
            if region not in self.structures_in_section:
                self.structures_in_section = np.hstack((self.structures_in_section,region))
        #print(structures_in_section)
        # Sort structure ids by ontology
        self.sorted_indx = np.argsort(np.array(list(map(lambda x:len(x),self.rsp.structure_tree.ancestor_ids(self.structures_in_section)))))

    def get_structure_ids(self, acrnms):
        regions = []
        for acrnm in acrnms:
            regions.append(self.acrnm_map[acrnm])
        return regions

    def getColorVolume(self, rgb_levels=DEFAULT_COLOR_VALUES):
        level_adjusted_arrays = []
        for colori, int_level in zip(['red', 'green', 'blue'], rgb_levels):
            colarray = np.clip(self.int_arrays[colori], a_min=int_level[0], a_max=int_level[1]) - int_level[0]
            colarray = (colarray * 255. / (int_level[1] - int_level[0])).astype('uint8')
            level_adjusted_arrays.append(colarray)
        return np.stack(level_adjusted_arrays, axis=-1)

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
        intensity_values_red = np.zeros((linepts.shape[0],160))
        intensity_values_green = np.zeros((linepts.shape[0],160))
        intensity_values_blue = np.zeros((linepts.shape[0],160))

        self.coords = {}
        self.linepts = linepts
        self.intensityValues = intensity_values_red

        for j in range(linepts.shape[0]):
            self.coords[j] = (linepts[j, 0], linepts[j, 1], linepts[j, 2])

            for k in range(-80,80):
                try:
                    intensity_values_red[j,k+80] = (self.volumeRed[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    intensity_values_green[j,k+80] = (self.volumeGreen[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])
                    intensity_values_blue[j,k+80] = (self.volumeBlue[int(linepts[j,0]),int(linepts[j,1]),int(linepts[j,2]+k)])

                except IndexError:
                    pass
        
        # display image
        #intensity_values *= 2.
        #self.trigScale = intensity_values.shape[0] / intensity_values.shape[1] 
        self.int_arrays = {}
        self.int_arrays['red'] = intensity_values_red
        self.int_arrays['green'] = intensity_values_green
        self.int_arrays['blue'] = intensity_values_blue
        self.volArray = self.getColorVolume()

if __name__ == '__main__':
    # command line inputs - mouse id
    args = parser.parse_args()
    #input_resampled = pathlib.Path(args.inputResampledImages)
    #output_dir_csv = pathlib.Path(args.annotationFileLocation)
    mouse_id = args.mouseID

    mask = volumeMask(mouse_id)