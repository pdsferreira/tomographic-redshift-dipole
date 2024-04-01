import numpy as np
import healpy as hp
from scipy.optimize import brute
from scipy.stats import norm
import matplotlib.lines as mlines
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from ligo.skymap.tool.ligo_skymap_contour import main
from ligo.skymap.postprocess.ellipse import find_ellipse
from ligo.skymap.kde import BoundedKDE
from astropy.io import fits
from scipy.interpolate import interp1d
import healpy as hp
import numpy as np
import json

def doppler_boost(beta_var, beta_lat, beta_long, vecs, z):
	'''
	beta_var = abolute value of the beta vector (v/c)
	beta_lat = latitude
	beta_long = longitude
	vecs = unitary direction vectors for each object
	z = z value for each object
	'''
	beta_vec = beta_var*hp.ang2vec(beta_long,beta_lat,lonlat=True).reshape(1,3)
	dotprocut = 1/(1+np.einsum('j,ij', beta_vec, vecs)) # equivalent to 1/(1+np.sum(vecs*beta_vec,axis=1))
	gamma = np.sqrt(1-(beta_var*beta_var)) 
	boostedvalues = gamma*dotprocut*z
	return boostedvalues

def doppler_deboost(beta_vec, vecs, z):
	'''
	beta_vec = beta vector (v/c)
	vecs = unitary direction vectors for each object
	z = z value for each object
	'''
	beta_var = np.linalg.norm(beta_vec)
	dotprocut = (1+np.einsum('j,ij', beta_vec, vecs))/np.sqrt(1-(beta_var*beta_var)) # equivalent to (1+np.sum(vecs*beta_vec,axis=1))
	boostedvalues = dotprocut*z
	return boostedvalues 	
		
def quantil_binning_number_list_s(catalog,n_bins):
	'''
	catalog = z distribution to be divided in bins
	n_bins = number of bins
	'''		
	bins = np.quantile(catalog,q=np.linspace(0,1,n_bins+1))
	bin_number_list = np.searchsorted(bins, catalog, side='left') # or np.digitize(z_array, bins, right=True)
	return np.int16(bin_number_list) # As the number of bins is small, we can use int16 precision to speedup the np.argwhere function that will be used latter in the estimator. As long as your number of bins is smaller than 32767 will should not change this line.
	
def z_dipole_estimator(z_array,z_hat,w_syst,n_bins,x_min=-0.0021,x_max=0.0021,y_min=-0.0021,y_max=0.0021,z_min=-0.0021,z_max=0.0021,iter1_step=0.00007,iter2_step=0.000035,iter3_step=0.00001133333,ns_2=4,ns_3=2,n_threads=1):
	'''
	This function can be used to estimate the dipole due to Doppler effect over the redshift distribution.
	The estimator consider a gradient descending approach over a chi squared grid with 3 iterations.
	For more details about the estimator check the paper at ArXiv XXXX.XXXX.
	
	z_array = array of z values. The array can have multiple sources (e.g. QSO NGC and QSO SGC), each one in a different row of the array.
	z_hat = unitary direction vectors for each object.
	n_bins = the number of bins to using during the estimator computation
	x,y,z_min and x,y,z_max = the ranges of the chi squared grid for the 1st iteration
	iter1_step, iter2_step and iter3_step = the step size of each iteration
	The range of the 2nd iteration is ns_2*iter1_step around the minimum of the 1st iteration. 
	The range of the 3rd iteration is ns_3*iter2_step around the minimum of the 2nd iteration.
	n_threads = number of threads used during minimization
	'''
	global chi2_binned
	n_sources = z_array.shape[0] # north and south data from sdss are considered as different sources as they have different z distributions (the monopole will be calculated separately)
	def chi2_binned(dipole):
		chi2_array = np.zeros((n_sources,n_bins)) # array to storage the chi squared of all bins for each source
		dipole_norm = np.linalg.norm(dipole) # norm of the dipole vector tested
		gamma = np.sqrt(1-(dipole_norm*dipole_norm))
		for s_idx in range(n_sources): # sources
			z_array_s = z_array[s_idx]
			z_hat_s = z_hat[s_idx]
			w_syst_s = w_syst[s_idx]
			values_deboosted = doppler_deboost(dipole, z_hat_s, z_array_s) # deboost of the z_array
			bin_number_list = quantil_binning_number_list_s(values_deboosted,n_bins) 
			gamma_times_dotprod_inv = (gamma/(1+np.einsum('j,ij', dipole, z_hat_s)))
			for idx in range(0, n_bins): # bins
				filter_ar = np.argwhere(bin_number_list == idx+1).T # boolean mask of the bin
				w_syst_bin = w_syst_s[filter_ar]
				avg_value = np.average(values_deboosted[filter_ar],weights=w_syst_bin) # deboosted monopole of the bin
				dipole_term = avg_value*gamma_times_dotprod_inv[filter_ar]
				w_sqrd = w_syst_bin**2
				diff = w_syst_bin*(z_array_s[filter_ar]-dipole_term) # squared difference
				chi2_array[s_idx][idx] = np.inner(diff,diff)/np.sum(w_sqrd) # weighted chi squared
		chi2 = np.sum(chi2_array) # sum of all chi squares, for all sources and bins
		return chi2
	# the next lines will run the temporary function "chi2_binned" over the grid of values for beta_x, beta_y and beta_z.
	rranges = (slice(x_min, x_max, iter1_step), slice(y_min, y_max, iter1_step),slice(z_min, z_max, iter1_step)) # 1st iteration
	res = brute(chi2_binned, rranges, full_output=True, finish=None, workers=n_threads) # minimizing the squared difference
	rranges = (slice(res[0][0]-(ns_2*iter1_step), res[0][0]+(ns_2*iter1_step), iter2_step), slice(res[0][1]-(ns_2*iter1_step), res[0][1]+(ns_2*iter1_step), iter2_step),slice(res[0][2]-(ns_2*iter1_step), res[0][2]+(ns_2*iter1_step), iter2_step)) # 2nd iteration
	res = brute(chi2_binned, rranges, full_output=True, finish=None, workers=n_threads) # minimizing the squared difference
	rranges = (slice(res[0][0]-(ns_3*iter2_step), res[0][0]+(ns_3*iter2_step), iter3_step), slice(res[0][1]-(ns_3*iter2_step), res[0][1]+(ns_3*iter2_step), iter3_step),slice(res[0][2]-(ns_3*iter2_step), res[0][2]+(ns_3*iter2_step), iter3_step)) # 3rd iteration
	res = brute(chi2_binned, rranges, full_output=True, finish=None, workers=n_threads) # minimizing the squared difference
	return res		
	

    
