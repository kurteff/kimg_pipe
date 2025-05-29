# Create a video of brain activity as
# gaussian blobs on mesh
# This one will use the high gamma data
# and the event file to find when to
# start.

## YOU NEED AN ENVIRONMENT WHERE BOTH MAYAVI AND MNE WORK 
# I recommend pip installing mne in your kimg_pipe environment,
# this is probably the easiest.
import mne
import matplotlib
from matplotlib import pyplot as plt
from pyface.api import GUI  # This is needed to fix a mayavi bug with the screenshot function
# see docs at https://github.com/enthought/mayavi/issues/711
from kimg_pipe import kimg_pipe
from kimg_pipe.plotting.ctmr_brain_plot import ctmr_gauss_plot, el_add
import os
import numpy as np
from mayavi import mlab
from matplotlib import cm

# For animations, from pycortex
linear = lambda x, y, m: (1.-m)*x + m*y
mixes = dict(
    linear=linear,
    smoothstep=(lambda x, y, m: linear(x,y,3*m**2 - 2*m**3)),
    smootherstep=(lambda x, y, m: linear(x, y, 6*m**5 - 15*m**4 + 10*m**3))
)

subj = 'TCH3'
hem = 'lh'
block = 'B1'
atlas = 'cvs_avg35_inMNI152'  # Not actually

imaging_dir = '/Users/liberty/Library/CloudStorage/Box-Box/NIH_ECoG_imaging/TCH_imaging'
ecog_dir = '/Users/liberty/Library/CloudStorage/Box-Box/NIH_ECoG/TCH_ECoG'
data_dir = f'{ecog_dir}/sub-{subj}/{subj}_{block}/'

# Initialize the patient
patient = kimg_pipe.freeCoG(subj=subj+'_complete', subj_dir=imaging_dir, hem=hem)

lh_pial = patient.get_surf(hem='lh')#, template=atlas)
rh_pial = patient.get_surf(hem='rh')#, template=atlas)
elecs = patient.get_elecs()['elecmatrix']
anat = patient.get_elecs()['anatomy']
elec_names = [a[0][0] for a in anat]

fif_name = 'ecog_hilbAA_70to150_8band_notch_car_log.fif'
cond = 'HilbAA_70to150_8band'

raw = mne.io.read_raw_fif(f'{data_dir}/{cond}/{fif_name}')
event = mne.read_events(f'{data_dir}/{subj}_{block}_sentence-eve.txt')

# Get the electrodes that match the raw data
elecs_plot = []
for ch in raw.info['ch_names']:
    elecs_plot.append(elecs[elec_names.index(ch),:])
elecs_plot = np.array(elecs_plot)

if hem == 'rh': # haven't tested
    start_stop_azimuth=(45, -45)
    start_stop_elevation=(60, 90)
elif hem == 'lh': # looks good, rotates from top oblique to back lateral, might be too much rotation
    start_stop_azimuth=(180-45, 180+45)
    start_stop_elevation=(60, 90)
else:  # stereo, haven't tested
    start_stop_azimuth=(90, 90)
    start_stop_elevation=(0, 0)
mix_type='smootherstep'
mix = mixes[mix_type]

epochs = mne.Epochs(raw, event, tmin=-0.5, tmax=2.0)
dat = epochs.get_data()
data = dat.mean(0)  # Take the average across repetitions of this event

vmin = -np.abs(data).max()*0.8
vmax = np.abs(data).max()*0.8
min_opacity = 0.25 # For the electrodes, minimum opacity for not very active elecs

# Get a mapper from data space to RGB for coloring the electrodes
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.RdBu_r)

nchans, ntimes = data.shape

# Loop through all the time points for this epoch
for mi, m in enumerate(np.linspace(0, 1, ntimes)):
    new_az = mix(start_stop_azimuth[0], start_stop_azimuth[1], m)
    new_el = mix(start_stop_elevation[0], start_stop_elevation[1], m)
    
    mlab.figure(fgcolor=(0, 0, 0), bgcolor=(0.5, 0.5, 0.5), size=(1200,900))
    # Plot the brains
    if (hem == 'lh') or (hem == 'stereo'):
        print('lh or stereo')
        mesh, mlab = ctmr_gauss_plot(lh_pial['tri'], lh_pial['vert'],
                                 elecs=elecs_plot,
                                 weights=data[:,mi], 
                                 new_fig=False,
                                 vmin=vmin,
                                 vmax=vmax,
                                 opacity=0.3, 
                                 )
    if (hem == 'rh') or (hem=='stereo'):
        print('rh or stereo')
        if hem == 'stereo':
            new_fig = False
        else:
            new_fig = False
        mesh, mlab = ctmr_gauss_plot(rh_pial['tri'], rh_pial['vert'],
                         elecs=elecs_plot,
                         weights=data[:,mi], 
                         new_fig=new_fig,
                         vmin=vmin,
                         vmax=vmax,
                         opacity=0.3, 
                         )

    # Plot the electrodes
    points = []
    for eidx, elec in enumerate(elecs_plot):
        if (hem=='lh') and (elec[0]>0):
            continue
        elif (hem=='rh') and (elec[0]<0):
            continue
        else:
            # Get the color for this data point, get rid of the alpha
            # because mayavi only takes RGB
            current_color = mapper.to_rgba(data[eidx,mi])[:3]
            #print(current_color)
            # Set the opacity to be between [min_opacity] and 1, based
            # on the data values for this time point. Larger positive
            # and negative values will be more opaque, values close
            # to zero will be more transparent
            op = np.min([min_opacity + np.abs(data[eidx,mi])/vmax, 1.0])
            points.append(el_add(np.atleast_2d(elec), color=current_color, 
                                 ambient=0.4225, specular = 0.333, 
                                 specular_power = 66, 
                                 diffuse = 0.6995,
                                 msize=3, opacity=op))

    print(new_az, new_el)
    GUI().process_events()
    mlab.view(new_az, new_el)

    # Save the frame as an image
    arr = mlab.screenshot(antialiased=True)
    plt.figure(figsize=(20, 10))
    plt.imshow(arr, aspect='equal')
    plt.axis('off')
    plt.show()
    if not os.path.isdir(f'{data_dir}/videos/'):
        os.mkdir(f'{data_dir}/videos')
    if not os.path.isdir(f'{data_dir}/videos/frames'):
        os.mkdir(f'{data_dir}/videos/frames')
    plt.savefig(f'{data_dir}/videos/frames/event_mean_all_{hem}_{mi:05d}.png')
    plt.close('all')
    mlab.close()


ffmpeg='/Users/liberty/opt/anaconda3/envs/mne/bin/ffmpeg'
movie_name = 'event_mean_32'
frame_rate = int(raw.info['sfreq']/4)
png_files = os.path.join(data_dir, 'videos', 'frames', f'event_mean_all_{hem}_%05d.png')
movie_file = os.path.join(data_dir, 'videos', movie_name+'.mp4')
os.system('%s -r %d -start_number 0 -i %s \
    -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 \
    -pix_fmt yuv420p %s'%(ffmpeg, frame_rate, png_files, movie_file))