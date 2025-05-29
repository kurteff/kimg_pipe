'''
utils.py (formerly img_pipe_tools.py or gkip.py)
Helper scripts for use with kimg_pipe
created by Lynn Kurteff for the Hamilton Lab
Updated 5/29/25
'''

from kimg_pipe import kimg_pipe
from kimg_pipe.plotting.ctmr_brain_plot import ctmr_gauss_plot
from kimg_pipe.plotting.ctmr_brain_plot import el_add
import nibabel
import numpy as np
import pandas as pd
import os
from PIL import ImageColor
import warnings
from pyface.api import GUI
import pylab as pl
from kimg_pipe.SupplementalFiles import FS_colorLUT
from kimg_pipe.img_pipe import remove_whitespace
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import time
import warnings

# # # # # # # # # # # # # # # # #
# DATA LOADING / PREPROCESSING  #
# # # # # # # # # # # # # # # # #
def nearest_electrode_vert(self, cortex_verts, elecmatrix, dist_thresh=4):
	''' Find the vertex on a mesh that is closest to the given electrode
	coordinates.
	
	Parameters
	----------
	cortex_verts : array-like
		[nvertices x 3] matrix of vertices on the cortical surface mesh
	elecmatrix : array-like
		[nchans x 3] matrix of 3D electrode coordinates 
	dist_thresh : in millimeters, the maximum distance we will allow for 
				  projection. This allows us to avoid projecting deep
				  white matter electrodes to the surface when they are too
				  far away

	Returns
	-------
	vert_inds : array-like
		Array of vertex indices that are closest to each of the 
		electrode 
	nearest_verts : array-like
		Coordinates for the nearest cortical vertices
	'''

	nchans = elecmatrix.shape[0]
	d = np.zeros((nchans, cortex_verts.shape[0]))

	# Find the distance between each electrode and all possible vertices
	# on the surface mesh
	for chan in np.arange(nchans):
		d[chan,:] = np.sqrt((elecmatrix[chan,0] - cortex_verts[:,0])**2 + \
							(elecmatrix[chan,1] - cortex_verts[:,1])**2 + \
							(elecmatrix[chan,2] - cortex_verts[:,2])**2)

	# Find the index of the vertex nearest to each electrode
	vert_inds = np.argmin(d, axis = 1)
	vert_dist = np.min(d, axis = 1)

	nearest_verts = cortex_verts[vert_inds,:]
	if dist_thresh is not None:
		print("Setting ")
		nbads = len(np.where(vert_dist > dist_thresh)[0])
		print(f"There are {nbads} electrodes with a distance greater than {dist_thresh} mm")
		nearest_verts[vert_dist > dist_thresh, :] = np.nan
		vert_inds[vert_dist > dist_thresh] = np.nan

	return vert_inds, nearest_verts
def calc_crtx_distance(patient, elec_xyz, hem='rh', clip_to='pial'):
	'''
	Calculates the distance of a given electrode from the cortical
	surface in mm, then returns this value.
	* clip_to: default to 'pial', surface you'd like to clip electrodes to. 'insula' currently only other
			   supported mode as the majority of insula elecs are >4mm from the pial surface.
	'''
	if clip_to == 'pial':
		cortex_verts = patient.get_surf(hem=hem)['vert']
	elif clip_to == 'insula':
		cortex_verts = patient.get_surf(hem=hem,roi="insula")['vert']
	if len(elec_xyz.shape) != 2:
		# expand dims on single elec
		elec_xyz = np.expand_dims(elec_xyz,axis=0)
	nchans = elec_xyz.shape[0]
	d = np.zeros((nchans, cortex_verts.shape[0]))
	for chan in np.arange(nchans):
		d[chan,:] = np.sqrt((
			elec_xyz[chan,0]-cortex_verts[:,0])**2+(elec_xyz[chan,1]-cortex_verts[:,1])**2+(elec_xyz[chan,2]-cortex_verts[:,2])**2
		)
	# vert_inds = np.argmin(d, axis=1)
	vert_dist = np.min(d, axis=1)
	# nearest_verts = cortex_verts[vert_inds,:]
	return vert_dist
def cat_anat(subjs, elecfile_prefix="TDT_elecs_all_warped",
	tch_imaging_path="only specify if including tch subjs"):
	'''
	Concatenates electrode anatomy across subjects.
	Only works in a warped space so it defaults to that.
	'''
	lh_anat, rh_anat = [],[]
	for s in subjs:
		if 'TCH' in s:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='lh',subj_dir=tch_imaging_path)
		else:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='lh')
		_,subj_anat = clip_hem_elecs(patient,hem='lh',elecfile_prefix=elecfile_prefix)
		if subj_anat.shape[0] != 0:
			lh_anat.append(subj_anat)
		if 'TCH' in s:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='rh',subj_dir=tch_imaging_path)
		else:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='rh')
		_,subj_anat = clip_hem_elecs(patient,hem='rh',elecfile_prefix=elecfile_prefix)
		if subj_anat.shape[0] != 0:
			rh_anat.append(subj_anat)
	all_lh_anat = np.vstack((lh_anat))
	all_rh_anat = np.vstack((rh_anat))
	return all_lh_anat,all_rh_anat
def cat_elecs(subjs,elecfile_prefix="TDT_elecs_all_warped",
	tch_imaging_path="only specify if including tch subjs"):
	'''
	Concatenates electrodes across subjects.
	Only works in a warped space so it defaults to that.
	'''
	lh_elecs, rh_elecs = [],[]
	for s in subjs:
		if 'TCH' in s:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='lh',subj_dir=tch_imaging_path)
		else:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='lh')
		subj_elecs,_ = clip_hem_elecs(patient,hem='lh',elecfile_prefix=elecfile_prefix)
		if subj_elecs.shape[0] != 0:
			lh_elecs.append(subj_elecs)
		if 'TCH' in s:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='rh',subj_dir=tch_imaging_path)
		else:
			patient = kimg_pipe.freeCoG(subj=f'{s}_complete',hem='rh')
		subj_elecs,_ = clip_hem_elecs(patient,hem='rh',elecfile_prefix=elecfile_prefix)
		if subj_elecs.shape[0] != 0:
			rh_elecs.append(subj_elecs)
	all_lh_elecs = np.vstack((lh_elecs))
	all_rh_elecs = np.vstack((rh_elecs))
	return all_lh_elecs,all_rh_elecs
def clip_4mm_elecs(patient, hem='rh', elecfile_prefix="TDT_elecs_all",
	elecmatrix=None, anatomy=None, debug=False, return_idxs=False, clip_to="pial"):
	'''
	Clips electrodes according to hemisphere. For subcortical/whitematter electrodes,
	prunes any that are >4mm from the surface of the cortex as these have a steep
	drop-off in high gamma activity: https://iopscience.iop.org/article/10.1088/1741-2560/13/5/056013/pdf
	* clip_to: default to 'pial', surface you'd like to clip electrodes to. 'insula' currently only other
			   supported mode as the majority of insula elecs are >4mm from the pial surface.
	'''
	if debug:
		print(patient.subj)
	if elecmatrix is None and anatomy is None:
		all_xyz = patient.get_elecs(elecfile_prefix=elecfile_prefix)['elecmatrix']
		all_anat = patient.get_elecs(elecfile_prefix=elecfile_prefix)['anatomy']
	else:
		all_xyz = elecmatrix
		all_anat = anatomy
	all_xyz, all_anat = clip_hem_elecs(patient, hem=hem, elecfile_prefix=elecfile_prefix,
		elecmatrix=all_xyz, anatomy=all_anat)
	elecmatrix, anatomy, idxs = [], [], []
	for i,a in enumerate(all_anat):
		roi = condense_roi(a[3][0])
		if roi in ['subcort','whitematter','outside_brain']:
			# Calculate distance to cortical surface
			distance = calc_crtx_distance(patient, all_xyz[i], clip_to=clip_to)[0]
			if debug:
				print(f"Electrode {a[0][0]} is {distance} mm from {clip_to}")
			if distance <= 4:
				elecmatrix.append(all_xyz[i])
				anatomy.append(a)
				idxs.append(i)
		else:
			elecmatrix.append(all_xyz[i])
			anatomy.append(a)
			idxs.append(i)
	if return_idxs:
		return np.array(elecmatrix), np.array(anatomy), idxs
	else:
		return np.array(elecmatrix), np.array(anatomy)
def clip_hem_elecs(patient,hem='rh',elecfile_prefix="TDT_elecs_all",
	elecmatrix=None,anatomy=None, return_idxs=False):
	'''
	Clips electrodes according to hemisphere. Returns elecmatrix and anatomy.
	'''
	if elecmatrix is None and anatomy is None:
		print("No elecmatrix/anat specified; clipping from all patient elecs...")
		all_xyz = patient.get_elecs(elecfile_prefix=elecfile_prefix)['elecmatrix']
		all_anat = patient.get_elecs(elecfile_prefix=elecfile_prefix)['anatomy']
	else:
		print("Clipping from specified subset of patient elecs...")
		all_xyz = elecmatrix
		all_anat = anatomy
	elecmatrix, anatomy, idxs = [], [], []
	for i,e in enumerate(all_xyz):
		if hem == 'rh':
			# Positive X only
			in_hem = e[0] >= 0
		else:
			# Negative X only
			in_hem = e[0] <= 0
		if in_hem:
			if hem == 'rh':
				if '_lh_' in all_anat[i][3][0] or 'left' in all_anat[i][3][0].lower():
					# Sometimes a positive X coordinate is still a LH elec because it crosses the midline
					# This and the following if statements relating to skip_append check
					# for that to prevent false positive errors.
					skip_append = True
				else:
					skip_append = False
			if hem == 'lh':
				if '_rh_' in all_anat[i][3][0] or 'right' in all_anat[i][3][0].lower():
					skip_append = True
				else:
					skip_append = False
			if not skip_append:
				elecmatrix.append(e)
				anatomy.append(all_anat[i])
				idxs.append(i)
	if return_idxs:
		return np.array(elecmatrix),np.array(anatomy), idxs
	else:
		return np.array(elecmatrix),np.array(anatomy)
def clip_roi_elecs(patient=None,elecmatrix=None,anatomy=None,include=[],exclude=[],
	hem='rh',elecfile_prefix="TDT_elecs_all",return_idxs=False):
	'''
	Uses freesurfer ROIs to return an elecmatrix and an anatomy for the specified freesurfer ROIs only
	* include, exclude: a list of freesurfer ROIs to include or exclude. Need to specify one or the other.
	See function condense_roi() for a list of possible rois.
	'''
	if type(include) == str or type(exclude) == str:
		warnings.warn("Please pass 'include' and 'exclude' as lists, not str. Attempting to automatically format as list now, but this may cause issues!")
		if len(include) > 0:
			include = list(include)
		if len(exclude) > 0:
			exclude = list(exclude)
	if len(include) > 0:
		inc = True
	elif len(exclude) > 0:
		inc = False
	else:
		raise Exception("Need to specify electrodes to include or exclude; neither were specified.")
	if patient != None:
		print("Loading all electrodes from patient")
		all_xyz, all_anat = clip_hem_elecs(patient,hem=hem,elecfile_prefix=elecfile_prefix)
	else:
		print("Using specified subset of patient electrodes")
		all_xyz = elecmatrix
		all_anat = anatomy
	elecmatrix, anatomy, idxs = [], [], []
	for i,a in enumerate(all_anat):
		roi = condense_roi(a[3][0])
		if inc:
			if roi in include:
				# print(patient.subj, roi)
				elecmatrix.append(all_xyz[i])
				anatomy.append(a)
				idxs.append(i)
		else:
			if roi not in exclude:
				elecmatrix.append(all_xyz[i])
				anatomy.append(a)
				idxs.append(i)
	if return_idxs:
		return np.array(elecmatrix),np.array(anatomy), idxs
	else:
		return np.array(elecmatrix),np.array(anatomy)
def clip_outside_brain_elecs(patient,elecmatrix=None,anatomy=None,
	hem='rh',force_include=[],elecfile_prefix="TDT_elecs_all",return_idxs=False):
	'''
	Cross-references patient electrodes with "IN_BOLT" textfile and returns an array
	excluding electrodes specified as outside the brain (i.e., excludes electrodes in
	that textfile). Defaults to checking all available electrodes for a patient but can
	be used with a subset of electrodes if `elecmatrix` and `anatomy` are passed as args.
	* force_include: a list of electrode names that are to be returned even if they are present
					 in the exclusionary "IN_BOLT" textfile.
	'''
	if type(force_include) == str:
		warnings.warn("Please pass force_include elecs as a list, not str. Attempting to automatically format as list now, but this may cause issues!")
		force_include = list(force_include)
	if elecmatrix is None and anatomy is None:
		print("Loading all electrodes from patient")
		all_xyz, all_anat = clip_hem_elecs(patient,hem=hem,elecfile_prefix=elecfile_prefix)
	else:
		print("Using a specified subset of patient electrodes")
		all_xyz = elecmatrix
		all_anat = anatomy
	ch_names = [a[0][0] for a in all_anat]
	subj_dir = patient.subj_dir
	subj = patient.subj
	in_bolt_fpath = os.path.join(subj_dir,subj,"elecs",f"{subj.replace('_complete','')}_IN_BOLT.txt")
	if not os.path.isfile(in_bolt_fpath):
		raise Exception(f"Could not locate {subj.replace('_complete','')}_IN_BOLT.txt. Please make sure this file exists before running this script.")
	elecs_in_bolt = np.loadtxt(in_bolt_fpath, dtype=str, skiprows=1)
	if len(elecs_in_bolt.shape) > 0:
		if elecs_in_bolt.shape[0] != 0:
			elecs_in_bolt = list(elecs_in_bolt)
			no_elecs_in_bolt = False
		else:
			no_elecs_in_bolt = True
	else:
		elecs_in_bolt = [str(elecs_in_bolt)]
		no_elecs_in_bolt = False
	if no_elecs_in_bolt:
		warnings.warn(f"Subject {subj.replace('_complete','')} has no electrodes in bolt, skipping...")
		if return_idxs:
			return all_xyz, all_anat, None
		else:
			return all_xyz, all_anat
	else:
		elecmatrix, anatomy, idxs, drop_elecs = [], [], [], []
		for i,ch in enumerate(ch_names):
			if ch in elecs_in_bolt and ch not in force_include:
				drop_elecs.append(ch)
			else:
				elecmatrix.append(all_xyz[i])
				anatomy.append(all_anat[i])
				idxs.append(i)
		if len(drop_elecs) > 0:
			print(f"Dropping {len(drop_elecs)} channels classified as outside the brain: {drop_elecs}...")
		if return_idxs:
			return np.array(elecmatrix), np.array(anatomy), idxs
		else:
			return np.array(elecmatrix), np.array(anatomy)
def color_by_roi(anat,mode='rgb_float',
	cmap = {
		'frontal' : '#4b72db', # blue
		'temporal' : '#84e16c', # green
		'parietal' : '#9c53b7', # purple
		'occipital' : '#d0c358', # yellow
		'precentral' : '#6ec2d5', # cyan
		'postcentral' : '#cb4779', # magenta
		'insula' : '#e14333', # red
		'subcort' : '#b0b0b0', # dark grey
		'whitematter' : '#5a5a5a', # light grey
		'outside_brain' : '#262626' # black
	}):
	'''
	Takes a list of electrodes and colors them by anatomy.
	* mode: choose from:
		- 'hex' : returns hex value
		- 'rgb_int' : returns integer based RGB value (0-255)
		- 'rgb_float' : returns float-based RGB value (0.-1.), the kimg_pipe el_add default
	'''
	colors = []
	for row in anat:
		ch_name = row[0]
		fs_roi = row[3][0]
		roi = condense_roi(fs_roi)
		if mode == 'hex':
			colors.append(cmap[roi])
		elif mode == 'rgb_int':
			colors.append(ImageColor.getcolor(cmap[roi],"RGB"))
		elif mode == 'rgb_float':
			rgb_tuple = ImageColor.getcolor(cmap[roi],"RGB")
			colors.append(tuple((rgb_tuple[0]/255,rgb_tuple[1]/255,rgb_tuple[2]/255)))
	return np.array(colors)
def condense_roi(fs_roi):
	'''
	Condenses freesurfer labels into a set of manageable ROIs for visualization/analysis:
	frontal, temporal, parietal, occipital, precentral, postcentral, insula, subcortical, and white matter.
	'''
	condensed_rois = {
		'frontal' : ['superiorfrontal','caudalmiddlefrontal','S_suborbital','S_orbital_lateral','S_orbital_med-olfact','S_orbital-H_Shaped','G_orbital',
					 'S_front_inf','S_front_middle','S_front_sup','G_front_inf-Opercular','G_front_inf-Orbital', 'parsopercularis','rostralmiddlefrontal',
					 'G_front_inf-Triangul','G_front_middle','G_front_sup','Lat_Fis-ant-Horizont','Lat_Fis-ant-Vertical', 'ctx_rh_S_precentral_inf-part',
					 'ctx_lh_G_rectus','ctx_rh_G_rectus','ctx_rh_G_and_S_transv_frontopol','ctx_lh_G_and_S_transv_frontopol', 'ctx_lh_S_precentral_inf-part',
					 'ctx_lh_superiorfrontal','ctx_rh_superiorfrontal'],
		'temporal' : ['S_temporal_inf','S_temporal_sup','S_temporal_transverse','S_collat_transv_ant',
					  'G_temp_sup-G_T_transv','G_temp_sup-Lateral','G_temp_sup-Plan_polar','G_temp_sup-Plan_tempo',
					  'G_temporal_middle','G_temporal_sup-Lateral','Pole_temporal','ctx_lh_G_temporal_inf',
					  'ctx_rh_G_temporal_inf','superiortemporal','middletemporal','bankssts','inferiortemporal'],
		'parietal' : ['superiorparietal','inferiorparietal','supramarginal','S_interm_prim-Jensen','G_pariet_inf-Angular','G_pariet_inf-Supramar','G_parietal_sup',
					  'ctx_lh_S_intrapariet_and_P_trans','ctx_rh_S_intrapariet_and_P_trans','ctx_lh_G_precuneus','ctx_rh_G_precuneus',
					  'ctx_lh_S_parieto_occipital','ctx_rh_S_parieto_occipital'],
		'occipital' : ['cuneus','S_oc_middle_and_Lunatus','S_oc-temp_med_and_Lingual','S_calcarine','G_oc-temp_lat-fusifor',
					   'G_oc-temp_med-Lingual','G_occipital_middle','Pole_occipital','ctx_lh_G_cuneus','ctx_rh_G_cuneus',
					   'ctx_lh_G_occipital_sup','ctx_rh_G_occipital_sup','ctx_lh_G_and_S_occipital_inf','ctx_rh_G_and_S_occipital_inf',
					   'ctx_lh_S_oc-temp_lat','ctx_rh_S_oc-temp_lat', 'ctx_lh_S_oc_sup_and_transversal', 'ctx_rh_S_oc_sup_and_transversal'],
		'precentral' : ['precentral','S_precentral-inf-part','S_central','G_precentral','G_and_S_subcentral','ctx_lh_S_precentral-sup-part',
						'ctx_rh_S_precentral-sup-part','precentral'],
		'postcentral' : ['postcentral','S_postcentral','G_postcentral','ctx_rh_G_and_S_paracentral','ctx_lh_G_and_S_paracentral','postcentral'],
		'insula' : ['S_circular_insula_ant','S_circular_insula_inf','S_circular_insula_sup',
					'G_Ins_lg_and_S_cent_ins','G_insular_short','Lat_Fis-post'],
		'subcort' : ['Right-Cerebellum-Cortex','Right-Hippocampus','Right-Amygdala','Left-Hippocampus','Left-Amygdala','Left-Caudate',
					 'G_and_S_cingul-Ant','G_and_S_cingul-Mid-Ant','G_oc-temp_med-Parahip','S_subparietal',
					 'ctx_lh_G_cingul-Post-ventral','ctx_rh_G_cingul-Post-ventral','ctx_lh_G_cingul-Post-dorsal',
					 'ctx_rh_G_cingul-Post-dorsal','Left-Thalamus','Right-Thalamus','ctx_lh_G_and_S_cingul-Mid-Post',
					 'ctx_rh_G_and_S_cingul-Mid-Post','Left-Putamen','Right-Putamen','ctx_rh_G_subcallosal','ctx_lh_G_subcallosal'],
		'whitematter' : ['ctx-rh-unknown','Unknown','WM-hypointensities','S_pericallosal','Right-Cerebral-White-Cortex',
						 'Left-Cerebral-White-Matter','Right-Cerebral-White-Matter','Left-Cerebral-White-Cortex',
						 'Left-Cerebral-Cortex','Right-Cerebral-Cortex','CC_Anterior'],
		'outside_brain' : ['Right-Inf-Lat-Vent','Left-Lateral-Ventricle','Right-Lateral-Ventricle','Left-Inf-Lat-Vent','Left-VentralDC','Right-VentralDC']
	}
	try:
		condensed_roi = [r for r in condensed_rois.keys() if fs_roi[7:] in condensed_rois[r]][0]
	except:
		try:
			condensed_roi = [r for r in condensed_rois.keys() if fs_roi in condensed_rois[r]][0]
		except:
			print(fs_roi)
			print(fs_roi[7:])
			raise Exception(f'ROI {fs_roi} missing from condensed_rois dict, please add it')
	return condensed_roi
def condense_roi_md(fs_roi):
	'''
	An alternate ROI compression scheme borrowed from Maansi Desai's code.
	I rewrote it to save space but I preserved all her
	freesurfer label-to-condensed label associations where possible.
	The one exception is the insula, where I split it up into more fine-grained labels.
	I also added more labels that are present in my dataset but not hers, trying as best
	as possible to preserve her initial categorization.
	'''
	condensed_rois = {
		# Subcort/WM
		'cerebellum' : ['Cerebellum-Cortex'],
		'amyg' : ['Amygdala'],
		'hippo' : ['G_oc-temp_med-Parahip','Hippocampus'],
		'thal' : ['Thalamus'],
		'cing' : ['G_subcallosal','G_and_S_cingul-Mid-Post','G_cingul-Post-dorsal','G_cingul-Post-ventral','G_and_S_cingul-Mid-Ant','G_and_S_cingul-Ant','S_subparietal'],
		'striatum' : ['Putamen','Caudate'],
		'wm' : ['Cerebral-White-Matter','Cerebral-Cortex','CC_Anterior','Cerebral-White-Cortex','S_pericallosal','WM-hypointensities'],
		'outside_brain' : ['Unknown','unknown','Inf-Lat-Vent','Lateral-Ventricle','VentralDC'],
		# Frontal
		'SFG' : ['superiorfrontal','G_front_sup','ctx_lh_superiorfrontal','ctx_rh_superiorfrontal'],
		'SFS' : ['S_front_sup'],
		'MFG' : ['front_middle','caudalmiddlefrontal','G_front_middle','rostralmiddlefrontal'],
		'MFS' : ['S_front_middle'],
		'IFG' : ['G_front_inf-Orbital','G_front_inf-Triangul','G_front_inf-Opercular','parsopercularis'],
		'IFS' : ['S_front_inf','S_precentral_inf-part'],
		'OFC' : ['S_suborbital','S_orbital_lateral','S_orbital_med-olfact','S_orbital-H_Shaped','G_orbital'],
		'front_operculum' : ['Lat_Fis-ant-Horizont','Lat_Fis-ant-Vertical'],
		'front_pole' : ['G_rectus','G_and_S_transv_frontopol'],
		# Precentral
		'preCG' : ['G_precentral','precentral'],
		'preCS' : ['S_precentral', 'S_precentral-sup-part','S_precentral-inf-part'],
		'CS' : ['S_central'],
		# Postcentral
		'postCG' : ['G_postcentral','postcentral'],
		'postCS' : ['S_postcentral'],
		# Paracentral
		'paracentral' : ['G_and_S_paracentral'],
		# Subcentral
		'subcentral' : ['G_and_S_subcentral'],
		# Temporal
		'STG' : ['G_temp_sup-Lateral','superiortemporal','G_temporal_sup-Lateral'],
		'STS' : ['S_temporal_sup','bankssts'],
		'MTG' : ['G_temporal_middle','middletemporal'],
		'ITG' : ['G_temporal_inf','inferiortemporal'],
		'ITS' : ['S_temporal_inf'],
		'PT' : ['G_temp_sup-Plan_tempo'],
		'PP' : ['G_temp_sup-Plan_polar'],
		'HG' : ['G_temporal_transverse','G_temp_sup-G_T_transv','S_temporal_transverse'],
		'temp_pole' : ['Pole_temporal'],
		# Parietal
		'supramar' : ['G_pariet_inf-Supramar','supramarginal'],
		'angular' : ['G_pariet_inf-Angular'],
		'SPL' : ['superiorparietal','G_parietal_sup','G_precuneus','S_parieto_occipital'],
		'IPL' : ['inferiorparietal'],
		'intraparietal' : ['S_interm_prim-Jensen','S_intrapariet_and_P_trans'],
		# Insula
		'insula_inf' : ['S_circular_insula_inf'],
		'insula_sup' : ['S_circular_insula_sup'],
		'insula_post' : ['G_Ins_lg_and_S_cent_ins','Lat_Fis-post'],
		'insula_ant' : ['S_circular_insula_ant','G_insular_short'],
		# Occipital (I'm not very good at occipital lobe anatomy so some of this may be off...)
		'cuneus' : ['cuneus','G_cuneus','S_oc_sup_and_transversal'],
		'lunate' : ['S_oc_middle_and_Lunatus'],
		'lingual' : ['S_oc-temp_med_and_Lingual','G_oc-temp_med-Lingual'],
		'calcarine' : ['S_calcarine'],
		'occ_temp' : ['G_oc-temp_lat-fusifor','S_oc-temp_lat','S_collat_transv_ant'],
		'occ_pole' : ['G_occipital_middle','Pole_occipital','G_occipital_sup','G_and_S_occipital_inf']
	}
	# Remove information that's redundant across hemispheres
	fs_roi_strip = fs_roi.replace("ctx_rh_","").replace("ctx_lh_","").replace("Left-","").replace("Right-","")
	try:
		condensed_roi = [r for r in condensed_rois.keys() if fs_roi_strip in condensed_rois[r]][0]
	except:
		print(fs_roi_strip)
		raise Exception(f'ROI {fs_roi_strip} missing from condensed_rois dict, please add it')
	return condensed_roi

def convert_elecs_to_inflated(patient,pial_elecs,hem='rh',use_pial=False, anat=None, warp=False):
	'''
	Takes electrodes from a regular pial space
	and maps them onto the inflated pial surface.
	'''
	if not warp:
		inf_pial = patient.get_surf(roi='inflated',hem=hem)
		pial_mesh = patient.get_surf(hem=hem)
	else:
		inf_pial, _ = load_template_brain(patient.subj_dir, inflated=True, hem=hem)
		pial_mesh,_ = load_template_brain(patient.subj_dir, inflated=False, hem=hem)
		patient = kimg_pipe.freeCoG("cvs_avg35_inMNI152", hem=hem)
	if use_pial:
		# Basically for legacy, don't use this like ever
		vert_inds, nearest_verts = patient.nearest_electrode_vert(pial_mesh['vert'],pial_elecs)
		inf_elecs = inf_pial['vert'][vert_inds,:]
	elif anat is not None:
		inf_elecs = np.zeros(pial_elecs.shape)
		for i,e in enumerate(pial_elecs):
			roi = anat[i][3][0].replace(f"ctx_{hem}_","")
			try:
				roi_mesh = patient.get_surf(roi=roi, hem=hem)
			except:
				# print(roi, hem)
				if "G_and_S" in roi:
					roi_mesh = patient.get_surf(roi=roi.replace("G_and_S","G&S"),hem=hem)
				elif "temporal" in roi:
					roi_mesh = patient.get_surf(roi=roi.replace("temporal","temp"),hem=hem)
				elif condense_roi(anat[i][3][0]) in ['subcort','whitematter','outside_brain']:
					print(
						f"Projecting electrode {anat[i][0][0]} to cortical surface (fs roi {anat[i][3][0]})"
						)
					roi_mesh = pial_mesh
				else:
					print(anat[i][3][0])
					raise Exception(f"Unsure how to rename {roi}")
			# These indicies are only appropriate for the sub-mesh, do not match pial/larger inflated meshes
			inds_roi, verts_roi = patient.nearest_electrode_vert(roi_mesh['vert'], np.expand_dims(e,axis=0))
			inds_pial, verts_pial =  patient.nearest_electrode_vert(pial_mesh['vert'], verts_roi)
			inf_elecs[i,:] = inf_pial['vert'][inds_pial,:]
	else:
		raise Exception(
			"For use_pial=False to work, you need to specify anat for the electrodes as well. You can get this using patient.get_elecs()['anatomy']."
			)
	return inf_elecs
def get_ch_names(patient=None,anat=None):
	'''
	Returns channel names from the electrode matrix 'anatomy'.
	'''
	if patient is not None:
		anat = patient.get_elecs()['anatomy']
	elif anat is None:
		raise Exception("Need to specify patient or anatomy")
	return [e[0] for e in anat[:,0]]
def get_rois(patient, hem="rh", elecfile_prefix="TDT_elecs_all"):
	'''
	Returns a list of possible ROIs for a patient
	'''
	elecs, anat = clip_hem_elecs(
		patient, hem=hem,
		elecmatrix=patient.get_elecs(elecfile_prefix=elecfile_prefix)['elecmatrix'],
		anatomy=patient.get_elecs(elecfile_prefix=elecfile_prefix)['anatomy']
	)
	elecs, anat = clip_roi_elecs(
		exclude=['whitematter','subcort','outside_brain'],
		elecmatrix=elecs, anatomy=anat, hem=hem
	)
	return list(np.unique([a[3][0].replace(f"ctx_{hem}_","") for a in anat]))
def load_curvature(patient,hem):
	'''
	Loads the white matter / grey matter curvature from the freesurfer data
	for plotting on an inflated cortex.
	Returns a weights array that you can add when calling ctmr_gauss_plot.
	Should have the same shape as your pial verts.
	'''
	curv_fpath = f'{patient.patient_dir}/surf/{hem}.curv'
	curv = nibabel.freesurfer.io.read_morph_data(curv_fpath)
	curv = curv/curv.max() # normalize
	return curv
def load_template_brain(imaging_path=None, template_name="cvs_avg35_inMNI152",inflated=False,hem='rh'):
	'''
	Loads the template brain for a given hemisphere.
	'''
	if imaging_path is not None:
		patient = kimg_pipe.freeCoG(subj=template_name, hem=hem, subj_dir=imaging_path)
	else:
		patient = kimg_pipe.freeCoG(subj=template_name, hem=hem)
	if inflated:
		pial = patient.get_surf(hem=hem,roi="inflated")
		curv = load_curvature(patient,hem)
		# threshold curvature to be binary
		curv[np.where(curv<=0)[0]] = -1
		curv[np.where(curv>0)[0]] = 1
		return pial, curv
	else:
		pial = patient.get_surf(hem=hem)
		return pial, None
def pop_device(patient,tgt_device,elecfile_prefix="TDT_elecs_all"):
	'''
	Returns a list of electrodes minus a certain device.
	Useful for patients with devices in both hemispheres.
	'''
	elecs = patient.get_elecs(elecfile_prefix=elecfile_prefix)['elecmatrix']
	anat = patient.get_elecs(elecfile_prefix=elecfile_prefix)['anatomy']
	ch_names = [e[0] for e in anat[:,0]]
	new_elecs = []
	new_anat = []
	pop_elecs = []
	pop_anat = []
	for i,c in enumerate(ch_names):
		if tgt_device not in c:
			new_elecs.append(elecs[i,:])
			new_anat.append(anat[i,:])
		else:
			pop_elecs.append(elecs[i,:])
			pop_anat.append(anat[i,:])
	return np.array(new_elecs), np.array(new_anat), np.array(pop_elecs), np.array(pop_anat)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PLOTTING FUNCTIONS 																		  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def plt_single_subject(subj,subj_dir,
	hem='rh', warp=False, inflated=False, clip_outside_brain=True, show_density=False, color_elecs=False,
	device_name=[], single_elec=[], roi = [], bgcolor=(1.,1.,1.), medial=False,
	labels=True, label_scale=4, msize=5, opacity=0.8, gsp=50, savefig=[]):
	'''
	Plots all electrodes (minus subcortical, whitematter, outside brain) for a single subject.
	This function is an amalgamation of many different plotting scripts that I combined on 8/9/23.
	Currently, show_density only works when inflated==False.
	* clip_oustide_brain: bool, defaults to True. If True, excludes electrodes marked as outside the brain in the "IN_BOLT" textfile.
	'''
	patient = kimg_pipe.freeCoG(subj=f"{subj}_complete", hem=hem, subj_dir=subj_dir)
	if warp:
		elecfile_prefix = "TDT_elecs_all_warped"
		if len(roi) > 0:
			distance=125
			mni_pt = kimg_pipe.freeCoG(subj="cvs_avg35_inMNI152", hem=hem)
			try:
				pial = mni_pt.get_surf(roi=roi, hem=hem)
			except:
				if "G_and_S" in roi:
					pial = mni_pt.get_surf(roi=roi.replace("G_and_S","G&S"), hem=hem)
				elif "temporal" in roi:
					pial = mni_pt.get_surf(roi=roi.replace("temporal","temp"), hem=hem)
				else:
					raise Exception(f"Unsure how to rename {roi}")
		else:
			distance=400
			pial, curv = load_template_brain(hem=hem,inflated=inflated)
	else:
		elecfile_prefix = "TDT_elecs_all"
		if len(roi) > 0:
			distance = 125
			try:
				pial = patient.get_surf(roi=roi, hem=hem)
			except:
				if "G_and_S" in roi:
					pial = patient.get_surf(roi=roi.replace("G_and_S","G&S"), hem=hem)
				elif "temporal" in roi:
					pial = patient.get_surf(roi=roi.replace("temporal","temp"), hem=hem)
				else:
					raise Exception(f"Unsure how to rename {roi}")
		else:
			distance = 400
			if inflated:
				pial = patient.get_surf(hem=hem, roi="inflated")
				curv = load_curvature(patient, hem)
				curv[np.where(curv<=0)[0]] = -1
				curv[np.where(curv>0)[0]] = 1
			else:
				pial = patient.get_surf(hem=hem)
	if len(device_name) > 0:
		_, _, elecs, anat = pop_device(patient, device_name, elecfile_prefix=elecfile_prefix)
	else:
		elecs, anat = clip_4mm_elecs(patient, hem=hem, elecfile_prefix=elecfile_prefix)
		if len(single_elec) > 0:
			if single_elec in [a[0][0] for a in anat]:
				elec_idx = [a[0][0] for a in anat].index(single_elec)
				elecs = np.expand_dims(elecs[elec_idx],axis=0)
				anat = np.expand_dims(anat[elec_idx],axis=0)
			else:
				raise Exception(
					f"Single electrode you're trying to plot ({single_elec}) was clipped because it is >4mm from cortical surface."
				)
	if len(elecs) > 0 and clip_outside_brain:
		elecs, anat = clip_outside_brain_elecs(patient,elecmatrix=elecs,anatomy=anat,hem=hem,elecfile_prefix=elecfile_prefix)
	if hem == 'rh':
		label_offset = 1.5
		if medial:
			azimuth = 175
		else:
			azimuth = -5
	else:
		label_offset = -1.5
		if medial:
			azimuth = -5
		else:
			azimuth = 175
	if len(elecs) > 0:
		if len(roi) > 0:
			roi_elecs, roi_anat = [], []
			for i,e in enumerate(elecs):
				fs_label = anat[i][3][0]
				if fs_label.replace(f"ctx_{hem}_","") == roi:
					roi_elecs.append(e)
					roi_elecs.append(anat[i])
			elecs = np.array(roi_elecs)
			anat = np.array(roi_anat)
			if warp:
				inds, verts = mni_pt.nearest_electrode_vert(pial['vert'], elecs)
			else:
				inds, verts = patient.nearest_electrode_vert(pial['vert'], elecs)
			elecs[:,] = pial['vert'][inds,:]
		if labels:
			labels = np.array([a[0][0] for a in anat])
		else:
			labels = None
		if color_elecs:
			elec_color = color_by_roi(anat)
		if inflated:
			elecs = convert_elecs_to_inflated(patient, elecs, hem=hem, anat=anat, warp=warp)
			mesh, mlab = ctmr_gauss_plot(tri=pial['tri'], vert=pial['vert'],
				brain_color=curv, cmap=cm.gray_r, vmin=-2, vmax=8, opacity=opacity, bgcolor=bgcolor)
		else:
			if show_density:
				weights = np.ones(elecs.shape[0])
				mesh, mlab = ctmr_gauss_plot(tri=pial['tri'], vert=pial['vert'], opacity=opacity,
					weights=weights, elecs=elecs, gsp=gsp, vmin=-1, vmax=1, cmap=cm.coolwarm, show_colorbar=False, bgcolor=bgcolor)
			else:
				mesh, mlab = ctmr_gauss_plot(tri=pial['tri'], vert=pial['vert'], opacity=opacity, bgcolor=bgcolor)
		if color_elecs:
			el_add(elecs, msize=msize, labels=labels, label_scale=label_scale, label_offset=label_offset, color=elec_color)
		else:
			el_add(elecs, msize=msize, labels=labels, label_scale=label_scale, label_offset=label_offset)
		if len(savefig) > 0:
			mlab.view(azimuth=azimuth, elevation=90, distance=distance)
			time.sleep(1)
			GUI().process_events()
			time.sleep(1)
			arr = mlab.screenshot(antialiased=True)
			fig = plt.figure(figsize=(20,20))
			arr, xoff, yoff = remove_whitespace(arr)
			pl.imshow(arr, aspect='equal')
			plt.axis('off')
			plt.tight_layout()
			plt.savefig(savefig, transparent=True)
			mlab.close()
			time.sleep(1)
			GUI().process_events()
			time.sleep(1)
			print("Screenshot saved to: ", savefig)
		return mlab
	else:
		return None
def plt_single_insula(subj,subj_dir,
	hem='rh', warp=False, clip_outside_brain=True, single_elec=[], show_density=False, color_elecs=False, bgcolor=(1.,1.,1.),
	medial=False,labels=True, label_scale=3, msize=3, opacity=0.8, gsp=50, savefig=[]):
	'''
	Yeah, the insula gets its own function :)
	Plots all electrodes in the insula on an insula mesh for a given subject.
	* clip_oustide_brain: bool, defaults to True. If True, excludes electrodes marked as outside the brain in the "IN_BOLT" textfile.
	'''
	patient = kimg_pipe.freeCoG(subj=f"{subj}_complete", hem=hem, subj_dir=subj_dir)
	if warp:
		elecfile_prefix = "TDT_elecs_all_warped"
		mni_pt = kimg_pipe.freeCoG(subj="cvs_avg35_inMNI152", hem=hem)
		pial = mni_pt.get_surf(hem=hem, roi='insula')
	else:
		elecfile_prefix = "TDT_elecs_all"
		pial = patient.get_surf(hem=hem, roi='insula')
	elecs, anat = clip_4mm_elecs(patient, hem=hem, elecfile_prefix=elecfile_prefix)
	elecs, anat = clip_roi_elecs(elecmatrix=elecs, anatomy=anat, include=['insula'])
	if len(single_elec) > 0:
		if single_elec in [a[0][0] for a in anat]:
			elec_idx = [a[0][0] for a in anat].index(single_elec)
			elecs = np.expand_dims(elecs[elec_idx],axis=0)
			anat = np.expand_dims(anat[elec_idx],axis=0)
		else:
			raise Exception(
				f"Single electrode you're trying to plot ({single_elec}) was clipped because it is >4mm from cortical surface."
			)
	if len(elecs) > 0 and clip_outside_brain:
		elecs, anat = clip_outside_brain_elecs(patient,elecmatrix=elecs,anatomy=anat,hem=hem,elecfile_prefix=elecfile_prefix)
	if hem == 'rh':
		label_offset = 1.5
		if medial:
			azimuth = 175
		else:
			azimuth = -5
	else:
		label_offset = -1.5
		if medial:
			azimuth = -5
		else:
			azimuth = 175
	if len(elecs) > 0:
		# Snap to nearest vertex on the insula mesh
		snap_elecs = np.zeros(elecs.shape)
		for i,e in enumerate(elecs):
			if warp:
				inds_roi, verts_roi = mni_pt.nearest_electrode_vert(pial['vert'],
					np.expand_dims(e,axis=0))
			else:
				inds_roi, verts_roi = patient.nearest_electrode_vert(pial['vert'],
					np.expand_dims(e,axis=0))
			snap_elecs[i,:] = pial['vert'][inds_roi,:]
		elecs = snap_elecs
		if labels:
			labels = np.array([a[0][0] for a in anat])
		else:
			labels = None
		if color_elecs:
			elec_color = color_by_roi(anat)
		if show_density:
			weights = np.ones(elecs.shape[0])
			mesh, mlab = ctmr_gauss_plot(tri=pial['tri'], vert=pial['vert'], opacity=opacity,
				weights=weights, elecs=elecs, gsp=gsp, vmin=-1, vmax=1, cmap=cm.coolwarm, show_colorbar=False, bgcolor=bgcolor)
		else:
			mesh, mlab = ctmr_gauss_plot(tri=pial['tri'], vert=pial['vert'], opacity=opacity, bgcolor=bgcolor)
		if color_elecs:
			el_add(elecs, msize=msize, labels=labels, label_scale=label_scale, label_offset=label_offset, color=elec_color)
		else:
			el_add(elecs, msize=msize, labels=labels, label_scale=label_scale, label_offset=label_offset)
		if len(savefig) > 0:
			mlab.view(azimuth=azimuth, elevation=90, distance=125)
			time.sleep(1)
			GUI().process_events()
			time.sleep(1)
			arr = mlab.screenshot(antialiased=True)
			fig = plt.figure(figsize=(20,20))
			arr, xoff, yoff = remove_whitespace(arr)
			pl.imshow(arr, aspect='equal')
			plt.axis('off')
			plt.tight_layout()
			plt.savefig(savefig, transparent=True)
			mlab.close()
			time.sleep(1)
			GUI().process_events()
			time.sleep(1)
			print("Screenshot saved to: ", savefig)
		return mlab
	else:
		return None