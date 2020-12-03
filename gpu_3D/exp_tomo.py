# General 3D - SIRT Reconstruction with Positivity Constraint. 

from pytvlib import *
import mpi_astra_ctvlib
import numpy as np
import time
from tqdm import tqdm
check_cuda()
########################################

# File Name
vol_size=''
file_name='bowtie.h5'

# Number of Iterations (Main Loop)
Niter = 200

alg = 'SIRT' # Algorithm
initAlg = '' # Algorithm Parameters (ie Projection Order or Filter)

# Descent Parameter and Reduction (ART)
beta0 = 0.5
beta_red = 0.995

# Save Final Reconstruction. 
saveRecon = True

##########################################

#Read Image. 
#(fName, tiltSeries) = load_data(vol_size,file_name)
#tiltAngles = np.load('Tilt_Series/'+fName+'_tiltAngles.npy')
(fName, tiltAngles, tiltSeries) = load_h5_data(vol_size,file_name)
(Nslice, Nray, Nproj) = tiltSeries.shape
Nproj = tiltAngles.shape[0]

# Initialize C++ Object.. 
tomo = astra_ctvlib.astra_ctvlib(Nslice, Nray, Nproj, np.deg2rad(tiltAngles))
initialize_algorithm(tomo, alg, initAlg)

# Create Projections Vector
b = np.zeros([tomo.NsliceLoc(), Nray*Nproj])
for s in range(tomo.NsliceLoc()):
    b[s,:] = tiltSeries[s+tomo.firstSlice(),:,:].transpose().ravel()
tomo.set_tilt_series(b)

dd_vec = np.zeros(Niter)
beta = beta0

if tomo.rank() == 0: print('Starting Reconstruction')

#Main Loop
for i in tqdm(range(Niter)): 

    run(tomo, alg, beta)

    # Data Distance
    dd_vec[i] = tomo.data_distance()

if tomo.rank() == 0:
	print('Reconstruction Complete, Saving Data..')
	print('Save Recon :: {}'.format(saveRecon))

#Save all the results to h5 file. 
fDir = 'results/' + fName + '_' + alg
meta = {'vol_size':vol_size,'Niter':Niter, 'initAlg':initAlg, 'beta':beta, 'beta_red': beta_red}
results = {'dd':dd_vec}
mpi_save_results([fDir, fName], tomo, saveRecon, meta, results)
