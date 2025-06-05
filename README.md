## `kimg_pipe` – kurteff's `img_pipe` ##
This is a fork of the Chang Lab's [img_pipe](https://github.com/ChangLabUcsf/img_pipe). Details about the original pipeline are published in [_Frontiers in Neuroinformatics_](https://doi.org/10.3389/fninf.2017.00062). The tool is great, so great in fact that I've been using it pretty consistently for pretty much my entire scientific career at this point. I've made some modifications to the source code in my local copy and written a bunch of helper scripts, and (WIP!!) more recently I have expanded the package so it can parse ROSA data in a fashion similar to the MATLAB-based [VERA](https://github.com/neurotechcenter/VERA/). I maintain this fork for my own uses and have no intent to publish any related content beyond what can be found in this repo, but if you find it useful, please let me know!

Maintained/forked by Lynn Kurteff, laboratories of Greg Hickok, UC Irvine & Liberty Hamilton, UT Austin

Email: lkurteff@uci.edu

Massive credit goes to Liberty Hamilton, an original `img_pipe` author and my PhD advisor. Even most of the "original" bits of code I've added in this fork she has had her hands in at some point. Other contributors from the Hamilton Lab to this repo in some form are Maansi Desai (PhD, now at Paradromics) and Alyssa Field (former Hamilton Lab RA). I'd also like to thank Markus Adamek, a contributor to VERA, for answering my questions via email and Mike Socha, a rep at Zimmer Bio, for answering my questions about how ROSA works via email and Zoom. I have tried to provide credit where deserved in the docstrings.

## Added features in `kimg_pipe`
* Handy plotting and analysis shortcuts
* Richer anatomical labeling
* Support for inflated meshes
* Standalone electrode picker (written by Liberty Hamilton)
* BIDS-compliant JSON generation (written by Liberty Hamilton)
* Animation rendering for presentations/supplements (written by Liberty Hamilton)

### `utils.py` – a one-stop shop for helper functions ###
This file contains a lot of helpful shortcuts for common analysis steps that I developed during my PhD work I will briefly describe the various top-down uses in this table:
<table>
    <tr>
        <th width=400px> Functionality </th>
        <th width=600px> Supporting Functions </th>
    </tr>
    <tr>
        <td>Snap electrode to the surface of a mesh</td>
        <td><pre lang="python">
kimg_pipe.utils.nearest_electrode_vert()</pre>
        </td>
    </tr>
    <tr>
        <td>Compute distance between electrode and cortical surface</td>
        <td><pre lang="python">
kimg_pipe.utils.calc_crtx_distance()</pre>
        </td>
    </tr>
    <tr>
        <td>Return electrode subsets based on common inclusion criteria</td>
        <td><pre lang="python">
kimg_pipe.utils.clip_hem_elecs()
kimg_pipe.utils.clip_4mm_elecs()
kimg_pipe.utils.clip_roi_elecs()
kimg_pipe.utils.clip_outside_brain_elecs()
kimg_pipe.utils.pop_device()</pre>
        </td>
    </tr>
    <tr>
        <td>Across-subject concatenation</td>
        <td><pre lang="python">
kimg_pipe.utils.cat_anat()
kimg_pipe.utils.cat_elecs()</pre>       
        </td>
    </tr>
    <tr>
        <td>Pick electrode colors for use with el_add</td>
        <td><pre lang="python">
kimg_pipe.utils.color_by_roi()</pre>
        </td>
    </tr>
    <tr>
        <td>Macro-anatomical labeling</td>
        <td><pre lang="python">
kimg_pipe.utils.condense_roi()
kimg_pipe.utils.condense_roi_md()</pre>
        </td>
    </tr>
    <tr>
        <td>Inflated surface generation and manipulation</td>
        <td><pre lang="python">
kimg_pipe.utils.convert_elecs_to_inflated()
kimg_pipe.utils.load_curvature()</pre>
        </td>
    </tr>
    <tr>
        <td>Shortcut parsers for freeCoG class and attributes</td>
        <td><pre lang="python">
kimg_pipe.utils.get_ch_names()
kimg_pipe.utils.get_rois()
kimg_pipe.utils.load_template_brain()</pre>
        </td>
    </tr>
    <tr>
        <td>Expanded vanilla plotting functions</td>
        <td><pre lang="python">
kimg_pipe.utils.plt_single_subject()
kimg_pipe.utils.plt_single_insula()</pre>
        </td>
    </tr>
</table> 

Oh yeah, I also expanded on the original workflow graphic to detail the changes in this fork:

![alt text](https://github.com/kurteff/kimg_pipe/blob/master/kimg_pipe/SupplementalFiles/kimg_workflow.png "kimg_pipe")



_(Yellow-bordered steps are still WIP)_

## Planned functionality for `kimg_pipe`
* Expanded `utils.py` functionality based on previously written but un-pushed spaghetti code _(partially implemented)_
* Shortcuts for mapping variables to 1D/2D colormaps, then to electrodes _(partially implemented)_
* Support for automatic electrode localization using Zimmer Bionet's ROSA surgical robot _(partially implemented)_
* Custom atlas support _(partially implemented)_
* Expanded BIDS utilities
* [MRIcroGL](https://www.nitrc.org/projects/mricrogl) implementation
* Better I/O with other commonly used Python neuroscience packages, such as the [NIPY](https://nipy.org/) package family and [ipyniivue](https://github.com/niivue/ipyniivue)

***

Below you will find the original readme that ships with img_pipe, which contains useful installation instructions among other handy resources:

# Original `README.md` #
## ![alt text](https://github.com/ChangLabUcsf/img_pipe/raw/master/img_pipe/SupplementalScripts/icons/leftbrain_blackbg.png "img_pipe") img_pipe: Image processing pipeline for ECoG data ![alt text](https://github.com/ChangLabUcsf/img_pipe/raw/master/img_pipe/SupplementalScripts/icons/rightbrain_blackbg.png "img_pipe") ##

[![Build Status](https://travis-ci.org/ChangLabUcsf/img_pipe.svg?branch=master)](https://travis-ci.org/ChangLabUcsf/img_pipe) [![PyPI version](https://badge.fury.io/py/img-pipe.svg)](https://badge.fury.io/py/img-pipe) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.996814.svg)](https://doi.org/10.5281/zenodo.996814)

Developed by Liberty Hamilton, David Chang, Morgan Lee at the Laboratory of Dr. Edward Chang, UC San Francisco
http://changlab.ucsf.edu

Email: edward.chang@ucsf.edu or liberty.hamilton@austin.utexas.edu with questions.

This contains the imaging pipeline as one importable python class for running a patient's
brain surface reconstruction and electrode localization/labeling.

The full capabilities of the pipeline are described in the paper: 
Hamilton LS, Chang DL, Lee MB, Chang EF (2017). [Semi-automated anatomical labeling and inter-subject warping of high-density intracranial recording electrodes in electrocorticography.](https://doi.org/10.3389/fninf.2017.00062) **Frontiers in Neuroinformatics** 11(62).

[Sample data is available on Zenodo](https://doi.org/10.5281/zenodo.996814), including an AC-PC aligned T1 scan, CT scan, and all intermediate and final files from the img_pipe processing pipeline.

## About ##
`img_pipe` is an open source python package for preprocessing of imaging data for use in intracranial electrocorticography (ECoG) and intracranial stereo-EEG analyses. This python package aims to provide a standardized interface for electrode localization, labeling, and warping to an atlas, as well as code to plot and display results on 3D cortical surface meshes. It gives the user an easy interface to create anatomically labeled electrodes that can also be warped to an atlas brain, starting with only a preoperative T1 MRI scan and a postoperative CT scan. 

Example results are shown below in the native subject space (left) and in the cvs_avg35_inMNI152 atlas space (right):

![alt text](https://github.com/ChangLabUcsf/img_pipe/raw/master/img_pipe/SupplementalFiles/img_pipe_results.png "img_pipe")

## Setup and Installation ##

To download this package, you will need:
* a MacOS or Linux machine (if you are using Windows, download a Linux Virtual Machine to use this package)
* __anaconda__ for Python version 2.7 or 3 (https://www.continuum.io/downloads)<br>
* __Freesurfer__ (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) version 5.3.0 or higher

After you download and install those dependencies, run the following commands in your terminal if using Python 2.7:

``` 
$ git clone https://github.com/changlabucsf/img_pipe
$ conda env create -f img_pipe/environment_py27.yml
$ source activate img_pipe_py2
$ ipython
$ import img_pipe
 ```

The following instructions should be used if you wish to work in Python 3:

```
$ git clone https://github.com/changlabucsf/img_pipe
$ conda env create -f img_pipe/environment_py35.yml
$ source activate img_pipe_py3
$ ipython
$ from img_pipe import img_pipe
```

After that, edit your ~/.bash_profile or ~/.bashrc and set the following environment variables with these lines:

```
export SUBJECTS_DIR=/path/to/freesurfer/subjects
export FREESURFER_HOME=/path/to/freesurfer/
source $FREESURFER_HOME/SetUpFreeSurfer.sh
```
Note that you can set `SUBJECTS_DIR` to wherever you want to place your subjects' imaging data - for example, `/Applications/freesurfer/subjects`.

Then in terminal, run `source ~/.bash_profile` or `source ~/.bashrc` to activate these environment variables.

To run `img_pipe`, you will need a high quality non-contrast T1 scan and a non-contrast CT scan. The T1 scan should ideally be 
AC-PC aligned before you start. Name the T1 scan T1.nii and place in `$SUBJECTS_DIR/your_subj/acpc`.  Name the CT scan CT.nii 
and place in `$SUBJECTS_DIR/your_subj/CT`.


You should now be able to import img_pipe from python. 
```python
>>> import img_pipe # Or in python 3, from img_pipe import img_pipe
>>> patient = img_pipe.freeCoG(subj='subject_name', hem='lh')
>>> patient.prep_recon()
>>> patient.get_recon()
```

If you have completed all of the steps, you can plot the brain with anatomically-labeled electrodes as follows:
```python
>>> import img_pipe
>>> patient = img_pipe.freeCoG(subj='subject_name', hem='lh')
>>> patient.plot_recon_anatomy()
```

Or just the brain with
```python
>>> patient.plot_brain()
```

The full workflow is shown as a flowchart below:

![alt text](https://github.com/ChangLabUcsf/img_pipe/raw/master/img_pipe/SupplementalFiles/workflow.png "img_pipe")

### Example images: ###
In addition to localization, labeling, and warping of electrodes, `img_pipe` includes some nice plotting functions for visualizing your data.  You can plot cortical and subcortical ROIs with different surface properties (opacity, color, specularity, etc) and you can plot electrodes either as spheres or as gaussian blobs on the cortical surface. Electrodes can be colored individually or all the same color. If you prefer to work in 2D (using matplotlib), there is a function `patient.auto_2D_brain()` that will create a brain image and corresponding projected 2D coordinates for even more flexibility (and some added speed since it avoids 3D rendering).

All of these images should be created after initializing the patient (as below).
```python
>>> import img_pipe
>>> patient = img_pipe.freeCoG(subj='subj_ID', hem='lh')
```

<table>
<tr>
<th width=300px>
Image
</th>
<th width=300px>
 Code
</th>
</tr>
<tr>
<td width=300px align="center">
 <img src="gallery/brain_left.png" alt="pial surface"/>
</td>
<td width=300px>
   <pre lang="python">
patient.plot_brain() 
   </pre>
</td>
</tr>
<tr>
<td align="center">
 <img src="gallery/brain_hipp.png" alt="pial surface with hippocampus"/>
</td>
<td>
   <pre lang="python">
pial = patient.roi('pial',opacity=0.3, representation='surface',
                   color=(0.9,0.8,0.5),gaussian=False)
hipp = patient.roi('lHipp', opacity = 0.8, representation = 'wireframe', 
                   color=(0.5, 0.3,0.5), gaussian=False) 
patient.plot_brain(rois=[pial, hipp])
   </pre>
</td>
</tr>
<tr>
<td align="center">
 <img src="gallery/ifg.png" alt="pial surface with IFG ROI"/>
</td>
<td>
   <pre lang="python">
roi_list = ['parsopercularis','parstriangularis','parsorbitalis']
patient.make_roi_mesh('pars',roi_list, showfig=False)
patient.plot_brain(rois=[patient.roi('pial',opacity=1.0),
                   patient.roi('lh_pars',color=(0.3,0.6,0.4))], 
                   screenshot=True, showfig=False)
   </pre>
</td>
</tr>
<tr>
<td align="center">
 <img src="gallery/recon_anatomy.png" alt="labeled electrodes on brain"/>
</td>
<td>
   <pre lang="python">
patient.plot_recon_anatomy()
   </pre>
</td>
</tr>
<tr>
<td align="center">
 <img src="gallery/brain_prepost.png" alt="electrodes on pre and postcentral gyri"/>
</td>
<td>
   <pre lang="python">
from plotting.ctmr_brain_plot import ctmr_gauss_plot
from plotting.ctmr_brain_plot import el_add
from matplotlib import cm
from matplotlib import pyplot as plt
patient = img_pipe.freeCoG(subj='EC108', hem='lh', subj_dir=subj_dir)
precentral_elecs = patient.get_elecs(roi='precentral')['elecmatrix']
postcentral_elecs = patient.get_elecs(roi='postcentral')['elecmatrix']
pial = patient.get_surf(hem='lh')
cmap = cm.Reds
precentral_colors = cmap(np.linspace(0,1,precentral_elecs.shape[0]))[:,:3]
cmap = cm.Blues
postcentral_colors = cmap(np.linspace(0,1,postcentral_elecs.shape[0]))[:,:3]
mesh, mlab = ctmr_gauss_plot(tri=pial['tri'], vert=pial['vert'])
el_add(precentral_elecs, color = precentral_colors)
el_add(postcentral_elecs, color = postcentral_colors)
   </pre>
</td>
</tr>
<tr>
<td align="center">
 <img src="gallery/atlaselecs.png" alt="warped electrodes on atlas brain"/>
</td>
<td>
   <pre lang="python">
subjs = ['EC108','EC125']
elecs = []

#get the electrode coordinate matrix for each subject
for s in subjs:
&nbsp;&nbsp;&nbsp;&nbsp;print s
&nbsp;&nbsp;&nbsp;&nbsp;patient = img_pipe.freeCoG(subj = s, hem = 'lh')
&nbsp;&nbsp;&nbsp;&nbsp;warped_elecs = patient.get_elecs(elecfile_prefix='TDT_elecs_all_warped')
&nbsp;&nbsp;&nbsp;&nbsp;elecs.append(warped_elecs['elecmatrix'])
#combine the electrode matrices from the different subjects into one matrix
elecmatrix = np.concatenate(elecs, axis=0)
#simply pass in the elecmatrix to plot_brain()
template = 'cvs_avg35_inMNI152'
atlas_patient = img_pipe.freeCoG(subj = template, hem='lh')
roi = atlas_patient.roi('pial', opacity=0.5)
atlas_patient.plot_brain(rois = [roi], 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;showfig=True, 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;screenshot=True, 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;elecs = elecmatrix,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;weights = None)
   </pre>
</td>
</tr>
</table>

If you find any bugs, please post in the Issues tab. 
