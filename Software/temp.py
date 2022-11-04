import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--inputResampledImages', help='Directory to resampeld images', required=True)
#parser.add_argument('-a', '--annotationFileLocation', help='Path for annotation csv file to be saved in this location', required=True)
parser.add_argument('--mouseID', help='Mouse ID of session')

if __name__ == '__main__':
    args = parser.parse_args()
    mouse_id = args.mouseID

    df_final_f1 = pd.read_csv('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe F1_channels_626791_warped.csv'.format(mouse_id))
    df_final_a1 = pd.read_csv('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/{}/Probe_A1_channels_626791_warped.csv'.format(mouse_id))
    reference = sitk.ReadImage('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/field_reference/average_template_25.nrrd')
    arr = sitk.GetArrayFromImage(reference).T

    print(df_final_f1)
    print(df_final_a1)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.scatter(df_final_f1['ML'] * 0.875, df_final_f1['DV'] * 3.84 * 0.875, s=15 ,alpha=0.95, c='red')
    ax.scatter(df_final_a1['ML'] * 0.875, df_final_a1['DV'] * 3.84 * 0.875, s=15 ,alpha=0.95)

    ax.imshow(arr.sum(axis=0), cmap='gray')

    plt.show()