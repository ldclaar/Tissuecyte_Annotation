# Tissuecyte_Annotation
App to annotate Tissuecyte Brains

# Environment Setup
First, create a new environment, and then activate the environment (app is currently using python 3.9):

```
conda create -n tissuecyte_annotation python=3.9.12
conda activate tissuecyte_annotation
```

After doing this, clone the repo:

```
git clone https://github.com/arjunsridhar12345/Tissuecyte_Annotation.git
```

And then `cd` to the directory where the repo was cloned and install the package and its dependencies

```
pip install -e .
```

# Using the App
### Getting the Affine Aligned Volume
First, the affine aligned volume needs to be extracted before proceeding to annotating. To do this, run the following command (example shown below)

```
python resample_images.py --user allen_institute_username --password allen_institute_password --mouseID mouseID
python resample_images.py --user arjun.sridhar --password arjun's password --mouseID 608671
```

Note: Access to HPC is a requirement as this script calls an executable on the network. To get access to HPC, go to ServiceNow and submit a ticket
This process can take some time (around 20-30 mins usually). The output will be saved to he following directory on the network

`/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID`

### Running the Annotation App
To use the annotation app, run the following command

```
python annotation_app_pyqtgraph_10.py --mouseID mouseID
```

The app may take several mintues to load, and once it has loaded, the following screen should be displayed

![image](https://github.com/arjunsridhar12345/Tissuecyte_Annotation/blob/main/images/annotation_app.png)
