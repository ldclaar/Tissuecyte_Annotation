# Tissuecyte_Annotation
App for Tissuecyte Annotation. Annotations are done in a 10 Micron space and then warped to the CCF space.

# Environment Setup
First, create a new environment, and then activate the environment (app is currently using python 3.9.12):

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
Then `cd Software` to begin using the app

# Using the App
### Getting the 10 Micron Volume
First, the 10 micron volume needs to be extracted before proceeding to annotating. To do this, run the following command giving the mouse id that needs to be annotated (example shown below)

```
python resample_images.py --user <allen_institute_username> --password <allen_institute_password> --mouseID <mouseID>
python resample_images.py --user arjun.sridhar --password arjun's password --mouseID 608671
```

Note: Access to HPC/Linux is a requirement as this script calls an executable on the network. To get access to HPC, go to ServiceNow and submit a ticket.
Running this process can take some time (around 20-30 mins usually). The output will be saved to he following directory on the network

`/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID`

### Running the Annotation App
To use the annotation app, run the following command

```
python annotation_app_pyqtgraph_10.py --mouseID <mouseID>
```

The app may take several mintues to load, and once it has loaded, the following screen should be displayed

![image](https://github.com/arjunsridhar12345/Tissuecyte_Annotation/blob/main/images/annotation_app.png)

### Functions of App
  * Point Lock - Disables annotating, deleting, and undoing points. When point lock is ON, all of these features are disabled
  * Annotating - Select the probe and day using the drop downs, and then disable point lock by clicking on the point lock button to set it to OFF. Then, click anywhre on the screen to create an annotation point. Once finished, enable point lock to prevent accidentally annotating unwanted points
  * Annotated points will be automatically saved to a CSV, which has 4 columns (AP, DV, ML, and probe_name). The CSV is stored in the following location
    `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/probe_annotations_mouseID.csv`
  * Undo Last Annotation - undos the last annotation made for the selected probe
  * Delete Points for Selected Probe - deletes all of the points for the selected probe on the current slice (only deletes all points on the current slice, not across all slices). A warning message will appear to confirm
  * Hide Points - Hides the points 
  * Show Points - Shows the points
  * RGB sliders/toggle - modifies the display for the given channels
  * Coronal, Horizontal, Sagittal buttons - change the view depending on which button is pressed
  * To move to different slices arcoss the brain, drag the main slider or use the left/right arrow keys
  * Zoom in/out - use mouse wheel
  * Can move image around my holding the left mouse button and then dragging anywhere on the screen

### Warp Points to CCF
  * Once all points have been annotated, click the Warp Annotations to CCF button to warp the points. Once this process has finished, a popup will appear saying that a new window will open with the warped points. Click OK and then wait for the CCF volume to be loaded. 
  * Then, scroll through the slices to see the warped points. Some of the warped points may appear on different slices, but they should generally be in the same regions of interest when compared with the 10 micron volume.
  * The warped points will be saved to following location
    `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/probe_annotations_mouseID_warped.csv`
  * In addition, each of the probes annotated will have a CSV file in the following location
    `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/Probe_probe_annotations_mouseID_warped.csv`
