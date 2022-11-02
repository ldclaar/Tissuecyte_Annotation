import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk

if __name__ == '__main__':
    df_final = pd.read_csv('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/608671/dummy_data/_alignment384.csv')
    reference = sitk.ReadImage('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/field_reference/average_template_25.nrrd')
    arr = sitk.GetArrayFromImage(reference).T

    print(df_final)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(df_final['ML'], df_final['DV'] * 3.84 * 0.875, s=15 ,alpha=0.95)

    ax.imshow(arr.sum(axis=0), cmap='gray')

    plt.show()