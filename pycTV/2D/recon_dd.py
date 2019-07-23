# General 2D - ART Reconstruction with Positivity Constraint.

import sys, os
sys.path.append('./Utils')
from matplotlib import pyplot as plt
from skimage import io
import numpy as np
import ctvlib 

algo = 'ART'
est = 'ones'

num_tilts = 30
beta_red = 0.995

if algo in 'ART':
    print('ART is Selected')
    Niter = 1000
    beta = 1.0
else:
    print('SIRT is Selected')
    Niter = 1000 
    beta = 0.0001

#Read Image. 
tiltSeries = io.imread('Test_Image/Co2P_256.tif')
tiltSeries = np.array(tiltSeries, dtype=np.float32)
tiltSeries /= np.amax(tiltSeries)
(Nx, Ny) = tiltSeries.shape
tiltSeries = tiltSeries.flatten()

# Generate Tilt Angles.
tiltAngles = np.linspace(0, 180, num_tilts, dtype=np.float32)

# Initialize C++ Object.. 
obj = ctvlib.ctvlib(Ny*num_tilts, Ny*Ny)

# Generate measurement matrix
A = obj.parallelRay(Ny, tiltAngles)
obj.rowInnerProduct()

b = np.transpose(A.dot(tiltSeries))
recon = np.ones([Nx, Ny], dtype=np.float32)
# recon = np.random.rand(Nx,Ny).astype(np.float32)
dd_vec = np.zeros(Niter)
recon_gif = np.zeros([Nx,Ny,Niter+1], dtype = np.float32)
recon_gif[:,:,0] = recon 

#Main Loop
for i in range(Niter): 

    if (i % 100 == 0):
        print('Iteration No.: ' + str(i+1) +'/'+str(Niter))
  
    if algo in 'ART':
        obj.ART(recon.ravel(), b, beta)
        beta *= beta_red
        # if i == Niter/2:
        # 	beta = 1.0
    else: 
        obj.SIRT(recon.ravel(), b, beta)
        beta *= beta_red

    #Positivity constraint 
    recon[recon < 0] = 0  

    recon_gif[:,:,i+1] = recon

    # DD Measurement
    g = A.dot(np.ravel(recon))
    dd_vec[i] = np.linalg.norm(g - b) / g.size

np.save(algo + '_' + est + '.npy', recon_gif)

x = np.arange(dd_vec.shape[0]) + 1

plt.figure(figsize=(5,4))
plt.plot(x,dd_vec,color='black', linewidth=2.0)
plt.title('DD', loc='left', fontweight='bold')
plt.title('Final DD: ' +str(dd_vec[-1]), loc='right')
plt.ylim(dd_vec[-1]/2,dd_vec[-1]*4)
plt.savefig(algo + '_' + est + '.png')
