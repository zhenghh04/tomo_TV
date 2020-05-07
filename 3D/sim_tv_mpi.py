# General 3D - ASD/TV Reconstruction with Positivity Constraint. 
# Intended for simulated datasets to measure RMSE and Volume's Original TV. 
# and to reconstruct large volume sizes (>1000^3) with Distributed Memory (OpenMPI)

import sys, os
from Utils.pytvlib import parallelRay, timer, load_data
#from mpi4py import MPI
import numpy as np
import Utils.mpi_ctvlib  as mpi_ctvlib
import time
########################################
import argparse
import h5py 

parser = argparse.ArgumentParser(description="Input for dynamic compressed sensing")
parser.add_argument("--niter", '-n', default=10, type=int, help='Number of iterations')
parser.add_argument("--data", '-d', default='Tilt_Series/256_au_sto.h5', type=str, help='Tilt data')
parser.add_argument("--measurement_matrix", '-m', default=None, type=str, help='Measurement Matrix')
parser.add_argument("--output", '-o', default="output.h5", help='output')
parser.add_argument("--noise", action='store_true', help='Whether to put noise or not')

args = parser.parse_args()

# Number of Iterations (Main Loop)
Niter = args.niter

# Number of Iterations (TV Loop)
ng = 10

# Parameter in ART Reconstruction.
beta = 0.25

# ART Reduction.
beta_red = 0.985

# Data Tolerance Parameter
eps = 0.019

# Reduction Criteria
r_max = 0.95
alpha_red = 0.95
alpha = 0.2

SNR = 100

#Outcomes:
noise = True                # Add noise to the reconstruction.
save_recon = 1           # Save final Reconstruction. 
##########################################

# Initalize pyMPI 
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()

#Read Image. (MPI_IO)
fv = h5py.File(args.data, 'r')
original_volume = fv['tiltSeries']
tiltAngles = fv['tiltAngles']
(Nslice, Nray, Nz) = original_volume.shape

# Generate Tilt Angles.
Nproj = tiltAngles.shape[0]
# Initialize C++ Object.. 
tomo_obj = mpi_ctvlib.mpi_ctvlib(Nslice, Nray, Nproj)
verbose = tomo_obj.get_rank()==0
if args.measurement_matrix==None:
    if tomo_obj.get_rank()==0:
        print("Generating measurement matrix")
    A = parallelRay(Nray, tiltAngles)
else:
    if verbose:
        print("reading measurement matrix from %s" %args.measurement_matrix)
    h5a = h5py.File(args.measurement_matrix, 'r')
    d = h5a["matrix"]
    A = np.zeros(d.shape)
    x, y = d.shape
    A = d[:x, :y]
    h5a.close()
if tomo_obj.get_rank()==0:
    print("Number of process: %d "%tomo_obj.get_nproc())
Nslice_loc = tomo_obj.get_Nslice_loc()
first_slice = tomo_obj.get_first_slice()

# Generate measurement matrix
if verbose:
    print("loading measurement matrix")
tomo_obj.load_A(A)
A = None
tomo_obj.rowInnerProduct()
if tomo_obj.get_rank()==0:
    print('Measurement Matrix is Constructed!')
# Load Volume and Collect Projections. 
t0 = time.time()
original_volume_loc  = np.zeros((Nslice_loc, Nray, Nz))
original_volume_loc = original_volume[first_slice:Nslice_loc+first_slice, :, :]
for s in range(Nslice_loc):
    tomo_obj.setOriginalVolume(original_volume_loc[s, :, :], s)
t1 = time.time()
if tomo_obj.get_rank()==0:
    print("loading local volume: ", t1 - t0)
# If creating simulation with noise, set background value to 1.
if noise:
    tomo_obj.set_background(1.0)
#   original_volume[original_volume == 0] = 1
if verbose:
    print("create projection")
tomo_obj.create_projections()

# Apply poisson noise to volume.
if noise:
    tomo_obj.poissonNoise(SNR)

#Measure Volume's Original TV
if verbose:
    print("compute original tv")
tv0 = tomo_obj.original_tv()

dd_vec, tv_vec = np.zeros(Niter), np.zeros(Niter)
rmse_vec, time_vec = np.zeros(Niter), np.zeros(Niter)

counter = 1 

t0 = time.time()

#Main Loop
if verbose:
    print("main loop")
from tqdm import tqdm
for i in tqdm(range(Niter)): 
    if ( i % 1 ==0 and tomo_obj.get_rank()==0):
        print('Iteration No.: ' + str(i+1) +'/'+str(Niter))
    if (verbose):
        print("iter %s"%i)
    tomo_obj.copy_recon()

    #ART Reconstruction. 
    tomo_obj.sART(beta, -1)

    #Positivity constraint
    tomo_obj.positivity()

    #ART-Beta Reduction
    beta *= beta_red 

    #Forward Projection.
    tomo_obj.forwardProjection(-1)

    #Measure Magnitude for TV - GD.
    if (i == 0):
        dPOCS = tomo_obj.matrix_2norm() * alpha
        dp = dPOCS / alpha
    else: # Measure change from ART.
        dp = tomo_obj.matrix_2norm() 

    # Measure difference between exp/sim projections.
    dd_vec[i] = tomo_obj.vector_2norm()

    #Measure TV.
    tv_vec[i] = tomo_obj.tv()

    #Measure RMSE.
    rmse_vec[i] = tomo_obj.rmse()

    tomo_obj.copy_recon() 

    #TV Minimization. 
    tomo_obj.tv_gd(ng, dPOCS)
    dg = tomo_obj.matrix_2norm()

    if(dg > dp * r_max and dd_vec[i] > eps):
        dPOCS *= alpha_red

    if (i+1)% 25 == 0:
        timer(t0, counter, Niter)

    counter += 1
    time_vec[i] = time.time() - t0

    if verbose:
        print("rmse: %s" %rmse_vec[i])

    #Save all the results to single matrix.
tomo_obj.save_recon(args.output, 0)
if tomo_obj.get_rank() == 0:
    f = h5py.File(args.output, 'w')
    results = np.array([dd_vec, eps, tv_vec, tv0, rmse_vec, time_vec])
#    os.makedirs('Results/'+ file_name +'_MPI/', exist_ok=True)
#    np.save('Results/' + file_name + '_MPI/results.npy', results)
    f['rmse'] = rmse_vec
    f['tv'] = tv_vec
    f['time'] = time_vec
    f['dd'] = dd_vec
    print("rmse_vec: ", rmse_vec)
    print("tv_vec: ", tv_vec)
    print("time_vec: ", time_vec)
    print("dd_vec: ", dd_vec)
    f.close()
#Get and save the final reconstruction.
#if parallel_hdf5 == False:
#    tomo_obj.gather_recon()
#    if save_recon and tomo_obj.get_rank()==0: 
#        recon = np.zeros([Nslice, Nray, Nray], dtype=np.float32, order='F')
#        #Use mpi4py to do MPI_Gatherv
#        for s in range(Nslice):
#            # recon_loc[s+first_slice,:,:] = tomo_obj.getRecon(s)
#            recon[s,:,:] = tomo_obj.getRecon(s)
#        f.create_dataset("recon", data=recon)
#        f.close()

tomo_obj.mpi_finalize()
