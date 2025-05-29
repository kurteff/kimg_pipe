from pyface.api import GUI  # This is needed to fix a mayavi bug with the screenshot function
# see docs at https://github.com/enthought/mayavi/issues/711

from kimg_pipe import kimg_pipe

subj='S0016'
hem='rh'
subj_dir='/Users/liberty/Library/CloudStorage/Box-Box/ECoG_imaging'

patient = kimg_pipe.freeCoG(subj=subj+'_complete', hem='stereo', subj_dir=subj_dir)

# First we want to get the mlab scene without having plotted it, otherwise
# when we call animate_scene it won't show up again for some reason
mlab = patient.plot_recon_anatomy(template='cvs_avg35_inMNI152', opacity=0.5, 
								  elecfile_prefix='TDT_elecs_all_warped',
								  showfig=False, screenshot=False)

# This will make a rotating brain animation based on the mlab scene above. If you
# wanted to animate and rotate a different scene, you'd just set that above instead.
patient.animate_scene(mlab[1], ffmpeg='/Users/liberty/opt/anaconda3/envs/mne/bin/ffmpeg', keep_frames=True,
					  nframes=200)

