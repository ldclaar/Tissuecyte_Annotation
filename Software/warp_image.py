import SimpleITK as sitk
import os
os.environ["OMP_NUM_THREADS"] = '1'
import sys
import subprocess
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import os
import argparse
import fnmatch

parser = argparse.ArgumentParser()
parser.add_argument('--volumeDir', help='Directory with saved volumes for annotated slices', required=True)
parser.add_argument('--annotationFile', help='CSV file with annotations', required=True)

# warps the volume using the reference, field and interpolator
def warp_execute( volume, reference, field, interpolator ) :
    deformWarped = sitk.Warp( volume, field, sitk.sitkLinear, 
                        reference.GetSize(), reference.GetOrigin(), reference.GetSpacing())
    return deformWarped

def get_points_line(nonzero_points, warped_dict, all_dict):
    for point in nonzero_points:
        ap = int(np.round(point[2]))
        dv = int(np.round(point[1]))
        ml = int(np.round(point[0]))

        warped_dict['AP'].append(ap)
        warped_dict['DV'].append(dv)
        warped_dict['ML'].append(ml)

        all_dict['AP'].append(ap)
        all_dict['DV'].append(dv)
        all_dict['ML'].append(ml)

# cluster based on number of annotations per probe
def cluster_annotations(k, nonzero_points, warped_dict, all_dict):
    print('K', k)
    kmeans = KMeans(n_clusters=k, random_state=0)
    labels = kmeans.fit_predict(nonzero_points)
    print(labels)
    centers = kmeans.cluster_centers_

    for center in centers:
        ap = int(np.round(center[2]))
        dv = int(np.round(center[1]))
        ml = int(np.round(center[0]))

        warped_dict['AP'].append(ap)
        warped_dict['DV'].append(dv)
        warped_dict['ML'].append(ml)

        all_dict['AP'].append(ap)
        all_dict['DV'].append(dv)
        all_dict['ML'].append(ml)
    """
    for i in range(k):
        indices = [j for j in range(len(labels)) if labels[j] == i] # get indices for each cluster
        points = []
          
        for index in indices: # get points for the cluster
            points.append(nonzero_points[index])

        points = np.array(points)
        #print(points)
        point = np.mean(points, axis=0) # take average

        ap = int(np.round(point[2]))
        dv = int(np.round(point[1]))
        ml = int(np.round(point[0]))

        warped_dict['AP'].append(ap)
        warped_dict['DV'].append(dv)
        warped_dict['ML'].append(ml)
    """

# helper function to warp the volume
def warp_volume_line(im, annotations, field, reference, probe, volume_dir, all_points, mouse_id):
    probe_file_name = probe.replace(' ', '_')
    warped_points = {'AP': [], 'DV': [], 'ML': []}
    #num_annotations = len(annotations) # get the number of points for the probe for the clustering
    

    volume_warped = warp_execute(im, reference, field, sitk.sitkLinear) # warp volume
    arr = sitk.GetArrayFromImage(volume_warped)
    nonzero_points = np.argwhere(arr > 0)
    num_annotations = len(nonzero_points)
    probe_column = [probe for i in range(num_annotations)]
    print(nonzero_points)

    #cluster_annotations(num_annotations, nonzero_points, warped_points, all_points)
    get_points_line(nonzero_points, warped_points, all_points)
    warped_points['probe_name'] = probe_column

    for probe in probe_column:
        all_points['probe_name'].append(probe)

    df_probe = pd.DataFrame(warped_points)
    df_probe.to_csv(os.path.join(volume_dir, '{}_annotations_{}_warped.csv'.format(probe_file_name, mouse_id)))

# helper function to warp the volume
def warp_volume(im, annotations, field, reference, probe, volume_dir, all_points, mouse_id):
    probe_file_name = probe.replace(' ', '_')
    warped_points = {'AP': [], 'DV': [], 'ML': []}
    num_annotations = len(annotations.loc[annotations['probe_name'] == probe].to_numpy()) # get the number of points for the probe for the clustering
    probe_column = [probe for i in range(num_annotations)]

    volume_warped = warp_execute(im, reference, field, sitk.sitkLinear) # warp volume
    arr = sitk.GetArrayFromImage(volume_warped)
    nonzero_points = np.argwhere(arr > 0)
    print(nonzero_points)

    cluster_annotations(num_annotations, nonzero_points, warped_points, all_points)
    warped_points['probe_name'] = probe_column

    for probe in probe_column:
        all_points['probe_name'].append(probe)

    df_probe = pd.DataFrame(warped_points)
    df_probe.to_csv(os.path.join(volume_dir, '{}_annotations_{}_warped_line.csv'.format(probe_file_name, mouse_id)))

# warps the points annotated based on the probe
def warp_points(output_dir, mouse_id, annotations, field, reference, progress):
    all_points = {'AP': [], 'DV': [], 'ML': [], 'probe_name': []}
    probes = annotations['probe_name'].unique()
    arr = np.zeros((528, 320, 456)) # 25 micron space
    
    steps = int(np.round(100 / len(probes))) # steps for progress bar
    t = 0
    for probe in probes:
        probe_file_name = probe.replace(' ', '_')
        print(probe)
        probe_annotations = annotations.loc[annotations['probe_name'] == probe]
        #print(probe_annotations)

        for annotation in probe_annotations.to_numpy():
            ap = int(np.round(annotation[0] / 2.5)) # affine aligned 25 micron space now
            dv = int(np.round(annotation[2] / 2.5)) # affine aligned 25 micron space now
            ml = int(np.round(annotation[1] / 2.5)) # affine aligned 25 micron space now
            arr[ap, dv, ml] = 1
        

        min_ind = int(np.round(probe_annotations['AP'].min() / 2.5))
        max_ind = int(np.round(probe_annotations['AP'].max() / 2.5))

        if min_ind == max_ind: # probes on same slice
            im = sitk.GetImageFromArray(arr[min_ind:min_ind+1, :, :].T)
        else: # probes across multiple slices
            im = sitk.GetImageFromArray(arr[min_ind:max_ind+1, :, :].T)

        im.SetSpacing((25, 25, 25))
        im.SetOrigin((min_ind * 25, 0, 0))
        warp_volume(im, annotations, field, reference, probe, output_dir, all_points, mouse_id)

        #sitk.WriteImage(im, os.path.join(volume_dir, probe_file_name + '.mhd'))
        print('Done')

        arr[arr > 0] = 0 # reset array to move on to next probe
        progress.setValue(t + steps) # update progress bar
        t += steps

    df_all_points = pd.DataFrame(all_points)
    df_all_points.to_csv(os.path.join(output_dir, 'probe_annotations_{}_warped.csv'.format(mouse_id)))
    progress.setValue(100)

def warp_lines(output_dir, mouse_id, linepts, field, reference):
    all_points = {'AP': [], 'DV': [], 'ML': [], 'probe_name': []}
    arr = np.zeros((528, 320, 456))
    for probe in linepts:
        points = linepts[probe]

        for annotation in points:
            ap = int(np.round(annotation[0])) # affine aligned 25 micron space now
            dv = int(np.round(annotation[1])) # affine aligned 25 micron space now
            ml = int(np.round(annotation[2])) # affine aligned 25 micron space now

            if ap < arr.shape[0] and dv < arr.shape[1] and ml < arr.shape[2]:
                arr[ap, dv, ml] = 1
        
        """
        min_ind = int(np.round(points[:, 0].min()))
        max_ind = int(np.round(points[:, 0].max()))

        if min_ind == max_ind: # probes on same slice
            im = sitk.GetImageFromArray(arr[min_ind:min_ind+1, :, :].T)
        else: # probes across multiple slices
            im = sitk.GetImageFromArray(arr[min_ind:max_ind+1, :, :].T)
        """
        im = sitk.GetImageFromArray(arr)
        im.SetSpacing((25, 25, 25))
        im.SetOrigin((0, 0, 0))
        warp_volume_line(im, points, field, reference, probe, output_dir, all_points, mouse_id)

    df_all_points = pd.DataFrame(all_points)
    df_all_points.to_csv(os.path.join(output_dir, 'probe_annotations_{}_warped_line.csv'.format(mouse_id)))

"""
def warp_volume(volume_dir, annotations, field, reference):
    volumes = os.listdir(volume_dir) # read in volumes
    volumes_mhd = fnmatch.filter(volumes, '*.mhd')
    #annotations = pd.read_csv(annotation_file)

    for volume_file in volumes_mhd:
        warped_points = {'AP': [], 'DV': [], 'ML': []}
        probe_file_name = volume_file[0:volume_file.index('.')]
        probe = probe_file_name.replace('_', ' ')
        print(probe)
        num_annotations = len(annotations.loc[annotations['probe_name'] == probe].to_numpy()) # get the number of probes for the clustering
        probe_column = [probe for i in range(num_annotations)]

        im = sitk.ReadImage(os.path.join(volume_dir, volume_file)) # read in volume for probe
        volume_warped = warp_slice(im, reference, field, sitk.sitkLinear) # warp volume
        arr = sitk.GetArrayFromImage(volume_warped)
        nonzero_points = np.argwhere(arr > 0)
        print(nonzero_points)

        cluster_annotations(num_annotations, nonzero_points, warped_points)
        warped_points['probe_name'] = probe_column

        df_probe = pd.DataFrame(warped_points)
        df_probe.to_csv(os.path.join(volume_dir, probe_file_name + '_annotations_warped.csv'))


    
    for volume in volumes_mhd:
        #vol_slice = volume[0: volume.index('_')]
        #view_slice = volume[volume.index('_'):2].upper()

        num_annotations = len(annotations.loc[annotations[view_slice] == vol_slice].to_numpy()) # get number of annotations in slice
        im = sitk.ReadImage(os.path.join(volume_dir, volume))

        volume_warped = warp_slice(im, reference, field, sitk.sitkLinear) # warp slice
        arr = sitk.GetArrayFromImage(volume_warped)

        nonzero_points = np.argwhere(arr > 0)

        # get warped points and add to dictionary
        cluster_annotations(num_annotations, nonzero_points, warped_points)

    annotations_sorted = annotations.sort_values('AP')
    warped_df = pd.DataFrame(warped_points)
    warped_df = warped_df.sort_values('AP')
    warped_df[['probe_name']] = annotations_sorted[['probe_name']].iloc[:, 0:4] # get probes for points

    warped_df.to_csv(os.path.join(volume_dir, 'warped_annotations.csv')) # save csv
    """
if __name__ == '__main__':
    args = parser.parse_args()
    vol_directory = args.volumeDir
    annotations = args.annotationFile

    #storage_directory = '//allen/scratch/aibstemp/arjun.sridhar/ForCorbett/image_series_1183416716'
    storage_directory = 'C:/608671/image_series_1183416716/local_alignment'
    #model_directory = '//allen/programs/celltypes/production/0378/informatics/model_september_2017/P56'
    model_directory = 'C:/608671/atlas'
    #input_file = os.path.join( working_directory, 'output_8', 'affine_vol_10.mhd')
    
    #affineAligned = sitk.ReadImage( input_file )
    #affineAligned.SetSpacing((25, 25, 25))
    field_file = os.path.join( storage_directory, 'dfmfld.mhd')
    reference_file = os.path.join( model_directory, 'average_template_25.nrrd')

    reference = sitk.ReadImage( reference_file )
    field = sitk.ReadImage( field_file )
    warp_volume(vol_directory, annotations, field, reference)
    #reference.SetOrigin((440, 0, 0))
    #field.SetOrigin((0, 0, 440))