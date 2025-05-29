# Plot sEEG electrodes on the MRI

import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('macosx')
import scipy.io
from matplotlib import cm
import os

class Patient:
    '''
    Patient class for seeg_mriplot. This allows you to get the MRI
    and electrode data for a given participant and plot electrodes
    in 2D sliceplanes with or without data overlaid.
    '''
    def __init__(self, subj, mri_dir, elecs_dir):
        '''
        Initialize the patient class
        Inputs:
            subj [str] : the subject ID 
            mri_dir [str] : path to the MRI directory containing brain.mgz
            elecs_dir [str] : path to the directory containing TDT_elecs_all.mat
        '''
        self.subj = subj.strip()  # remove leading and trailing spaces
        self.mri_dir = mri_dir.strip()  # remove leading and trailing spaces
        self.elecs_dir = elecs_dir.strip()  # remove leading and trailing spaces
        self.get_devices()
        self.load_mri()

    def get_devices(self):
        '''
        Get a list of the possible electrode devices and the number
        of electrodes in each. This loops through the electrode montage
        file to find this information.
        '''
        print("Getting a list of electrode devices and number of electrodes in each...")
        self.load_elecs()
        all_devices, inds, counts = np.unique([m.strip('1234567890') for m in self.montage], 
                                        return_counts=True, return_index=True)

        self.devices = dict()
        for i, d in enumerate(all_devices):
            self.devices[d] = [inds[i], counts[i]]  # Start index and # elecs

    
    def get_elecs(self, device_list):
        '''
        Get the electrode numbers for a given set of devices.
        '''
        elec_nums = []
        for device in device_list:
            start_elec = self.devices[device][0]
            num_elecs = self.devices[device][1]
            elec_nums += list(np.arange(start_elec, start_elec+num_elecs))
        return elec_nums

    def plot_close_electrodes(self):
        '''
        Plot pairs of electrodes that share at least one center 
        of mass plane of section
        '''
        com = dict() # Center of mass for each device
        for d in self.devices.keys():
            start_elec = self.devices[d][0]
            num_elecs = self.devices[d][1]
            com[d] = self.elecmatrix_RAS[start_elec:start_elec+num_elecs,:].mean(0)
        self.com = com

        #com_array = np.vstack((patient.com.values()))
        all_devices = list(self.devices.keys())

        close_pairs = []
        for i1, d1 in enumerate(all_devices):
            for i2, d2 in enumerate(all_devices):
                if i1 < i2:
                    if np.sum(np.abs(com[d1] - com[d2]) < 2) >= 2:
                        print(d1, d2)
                        close_pairs.append([d1, d2])

        for device_list in close_pairs:
            print(device_list)
            elec_nums = self.get_elecs(device_list)
            self.plot_device(device_name=None, elec_nums=elec_nums)


    def plot_all(self):
        '''
        Plot all of the electrode devices, each as a separate figure.
        '''
        for d in self.devices.keys():
            self.plot_device(d)


    def plot_device(self, device_name, elec_nums=None, data=None, vmin=None, vmax=None, 
                    cmap=cm.RdBu_r, crop=35, savefig=True):
        '''
        Plot an electrode device or a list of electrodes. If plotting a list,
        it's best if these match location in at least one plane, otherwise
        the plot is not interpretable.

        Input:
            device_name [str] : name of the device to plot, which can be found from
                                patient.get_devices() as the key values in the dict
                                For example, 'AIF-OF'
            elec_nums [list] :  List of electrode numbers, if you prefer to use these
                                instead of the device_name. If device_name is provided,
                                this elec_nums list is ignored. If you wish to use elec_nums,
                                you must set device_name=None
            data [np.array] :   1D vector of values for plotting on the electrodes. The length
                                of the vector should match the number of electrodes being
                                plotted.
            vmin [float] :      The minimum for the colorbar. If none, defaults to -abs(data).max()
            vmax [float] :      The maximum for the colorbar. If none, defaults to abs(data).max()
            cmap [matplotlib] : Matplotlib colormap. For example cm.RdBu_r or cm.Reds
            crop [int] :        Number of slices to crop around the MRI to decrease the amount
                                of empty space. Default = 35, you may need to play around.
            savefig [bool] :    Save figure (True/False).
        '''

        if device_name is not None:
            if device_name=='all':
                elecs = self.elecmatrix_RAS
                anat = self.anatomy
                msg = '\nPLOTTING ALL NOT RECOMMENDED - INTERPRET AS GLASS BRAIN'
            else:
                ind = self.devices[device_name][0]
                num_elecs = self.devices[device_name][1]
                msg = ''

                # Get the electrodes to plot
                elecs = self.elecmatrix_RAS[ind:ind+num_elecs,:]

                # Get the anatomical labels
                anat = self.anatomy[ind:ind+num_elecs]
        else:
            print("Overriding device and using electrode numbers for %d elecs"%(len(elec_nums)))
            elecs = self.elecmatrix_RAS[elec_nums, :]
            anat = [self.anatomy[e] for e in elec_nums]
            print(elecs.shape)

        elecs = elecs - crop

        subj_dat = self.mri.get_data()
        subj_dat = subj_dat[crop:-crop, crop:-crop, crop:-crop]

        # Find center of mass x of the electrode device so we can choose
        # MRI slices that best show the full device at once. Note that this
        # will not be the "true" location of the electrodes, it's simply
        # a representation!
        x = int(elecs[:,0].mean())
        y = int(elecs[:,1].mean())
        z = int(elecs[:,2].mean())

        # Get colors from the data, if provided
        if data is not None:
            data = data.ravel()
            if vmin is None:
                vmin=-np.abs(data).max()
            if vmax is None:
                vmax=np.abs(data).max()
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

            data_colors = []
            for d in data:
                data_colors.append(mapper.to_rgba(d))

        # Plot the brain and the electrodes on it
        ax1=plt.subplot(1,3,1).axes
        plt.imshow(subj_dat[x,:,:], cmap=cm.gray, vmin=0, vmax=subj_dat.max())
        if data is None:
            plt.plot(elecs[:,2], elecs[:,1], '.', markersize=6, color='r')
        else:
            for e in np.arange(elecs.shape[0]):
                print(e)
                plt.plot(elecs[e,2], elecs[e,1], '.', markersize=6, color=data_colors[e])
        plt.axis('equal'); ax = plt.gca(); ax.set_axis_off()

        ax2=plt.subplot(1,3,2).axes
        plt.imshow(subj_dat[:,y,:].T, cmap=cm.gray, vmin=0, vmax=subj_dat.max())
        if data is None:
            plt.plot(elecs[:,0], elecs[:,2], '.', markersize=6, color='r')
        else:
            for e in np.arange(elecs.shape[0]):
                plt.plot(elecs[e,0], elecs[e,2], '.', markersize=6, color=data_colors[e])
        plt.axis('equal'); ax = plt.gca(); ax.set_axis_off()
        
        ax3=plt.subplot(1,3,3).axes
        plt.imshow(subj_dat[:,:,z].T, cmap=cm.gray, vmin=0, vmax=subj_dat.max())   
        if data is None:
            plt.plot(elecs[:,0], elecs[:,1], '.', markersize=6, color='r')
        else:
            for e in np.arange(elecs.shape[0]):
                plt.plot(elecs[e,0], elecs[e,1], '.', markersize=6, color=data_colors[e])
        plt.axis('equal'); ax = plt.gca(); ax.set_axis_off()

        if device_name is not None:
            plt.suptitle('Device %s%s'%(device_name, msg))

        if savefig:
            plt.savefig(os.path.join(self.elecs_dir, '%s_%s_elecs.pdf'%(self.subj, device_name)))
        plt.show()


    def load_mri(self, mri_file='brain.mgz'):
        '''
        Load the MRI image from the MRI directory. Default is 'brain.mgz'
        but could also be 'aseg.mgz' or other images.
        '''
        print("Loading MRI image")
        self.mri = nib.freesurfer.load(os.path.join(self.mri_dir, mri_file))


    def load_elecs(self, elecs_file='TDT_elecs_all.mat'):
        '''
        Load electrodes from [elecs_file]
        Inputs:
            elecs_file [str] : Name of your electrodes file. This must be a mat file
                               with a variable 'anatomy' and a variable 'elecmatrix'.
                               This should have been created by kimg_pipe
        '''
        elecdata = scipy.io.loadmat(os.path.join(self.elecs_dir, elecs_file))
        self.anatomy = [elecdata['anatomy'][a][3][0] for a in np.arange(elecdata['anatomy'].shape[0])]
        self.montage = [elecdata['anatomy'][a][0][0] for a in np.arange(elecdata['anatomy'].shape[0])]
        self.elecmatrix = elecdata['elecmatrix']

        # Convert the electrode matrix into RAS coordinates, which is the
        # space needed to plot them on the MRI scan instead of on the surface.

        fsVox2RAS = np.array( [[-1., 0., 0., 128.], 
                               [0., 0., 1., -128.],
                               [0., -1., 0., 128.], 
                               [0., 0., 0., 1.]])


        # Convert surface RAS to voxel CRS
        elec = np.hstack((self.elecmatrix, np.ones((self.elecmatrix.shape[0], 1)) ))

        VoxCRS = np.dot(np.linalg.inv(fsVox2RAS), elec.transpose()).transpose()
        elec = VoxCRS
        elec = elec[:,0:3]

        self.elecmatrix_RAS = elec

if __name__ == "__main__":
    subj = input("Enter subject ID (e.g. S0003): ")
    mri_dir = input("Enter path to MRI directory (e.g. /Applications/freesurfer/subjects/S0003/mri )")
    elecs_dir = input("Enter path to elecs directory (e.g. /Applications/freesurfer/subjects/S0003/elecs )") 
    patient = Patient(subj, mri_dir, elecs_dir)
    patient.plot_all()

