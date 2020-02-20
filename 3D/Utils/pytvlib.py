from skimage import io
import scipy.ndimage
import numpy as np
import time
from tqdm import tqdm
import h5py

def tv_derivative(recon):
    r = np.lib.pad(recon, ((1, 1), (1, 1), (1, 1)), 'edge')
    v = (3 * r - np.roll(r, 1, axis=0) - \
                              np.roll(r, 1, axis=1) - np.roll(r, 1, axis=2)) / \
        np.sqrt(1e-8 + (r - np.roll(r, 1, axis=0))**2 + (r -
                  np.roll(r, 1, axis=1))**2 + (r - np.roll(r, 1, axis=2))**2)
    v += (r - np.roll(r, -1, axis=0)) / \
         np.sqrt( 1e-8 + ( np.roll(r, -1, axis=0) - r )**2 +
            ( np.roll(r, -1, axis=0) - np.roll(r, [-1,1], axis=(0,1)) )**2 +
            ( np.roll(r, -1, axis=0) - np.roll(r, [-1,1], axis=(0,2)) )**2 )
    v += (r - np.roll(r, -1, axis=1)) / \
         np.sqrt( 1e-8 + 
                 (np.roll(r, -1, axis=1) - np.roll(r, [-1,1], axis=(1,0)))**2 +
                 (np.roll(r, -1, axis=1) - r)**2 + 
                 ( np.roll(r, -1, axis=1) - np.roll(r, [-1,1], axis=(1,2)))**2) 
    v += (r - np.roll(r, -1, axis=2)) / \
         np.sqrt(1e-8 + 
                 (np.roll(r, -1, axis=2) - np.roll(r, [-1,1], axis=(2,0)))**2 +
                 (np.roll(r, -1, axis=2) - np.roll(r, [-1,1], axis=(2,1)))**2 +
                 (np.roll(r, -1, axis=2) - r)**2)
    v = v[1:-1, 1:-1, 1:-1]
    return v


def tv(recon):
    r = np.lib.pad(recon, ((1, 1), (1, 1), (1, 1)), 'edge')
    tv = np.sqrt(1e-8 + (r - np.roll(r, 1, axis=0))**2 +
                 (r - np.roll(r, 1, axis=1))**2 +
                 (r - np.roll(r, 1, axis=2))**2)
    tv = np.sum(tv[1:-1, 1:-1, 1:-1])
    return tv

def timer(t0, counter, Niter):
    timeLeft = (time.time() - t0)/counter * (Niter - counter)
    timeLeftMin, timeLeftSec = divmod(timeLeft, 60)
    timeLeftHour, timeLeftMin = divmod(timeLeftMin, 60)
    print('Estimated time to complete: %02d:%02d:%02d' % (timeLeftHour, timeLeftMin, timeLeftSec))

def generate_tilt_series(volume, angles, num_tilts):

    Ny = volume.shape[1]
    Nz = volume.shape[2]
    N = volume.shape[1]

    # pad volume
    pad_y_pre = int(np.ceil((N - Ny) / 2.0))
    pad_y_post = int(np.floor((N - Ny) / 2.0))
    pad_z_pre = int(np.ceil((N - Nz) / 2.0))
    pad_z_post = int(np.floor((N - Nz) / 2.0))
    volume_pad = np.lib.pad(
        volume, ((0, 0), (pad_y_pre, pad_y_post), (pad_z_pre, pad_z_post)),
        'constant')

    Nslice = volume.shape[0]  # Number of slices along rotation axis.
    tiltSeries = np.empty([Nslice, N, num_tilts], dtype=float, order='F')

    for i in range(num_tilts):
        # Rotate volume about x-axis
        rotatedVolume = np.empty_like(volume_pad)
        scipy.ndimage.interpolation.rotate(
            volume_pad, angles[i], axes=(1, 2), reshape=False, order=1,
            output=rotatedVolume)
        # Calculate projection
        tiltSeries[:, :, i] = np.sum(rotatedVolume, axis=2)
    return tiltSeries

def parallelRay(Nside, angles):

    Nray = Nside

    pixelWidth = 1.0
    rayWidth = 1.0

    # Suppress warning messages that pops up when dividing zeros
    np.seterr(all='ignore')
    Nproj = angles.size # Number of projections

    # Ray coordinates at 0 degrees.
    offsets = np.linspace(-(Nray * 1.0 - 1) / 2,
                          (Nray * 1.0 - 1) / 2, Nray) * rayWidth
    # Intersection lines/grid Coordinates
    xgrid = np.linspace(-Nside * 0.5, Nside * 0.5, Nside + 1) * pixelWidth
    ygrid = np.linspace(-Nside * 0.5, Nside * 0.5, Nside + 1) * pixelWidth

    # Initialize vectors that contain matrix elements and corresponding
    # row/column numbers
    rows = np.zeros((2 * Nside * Nproj * Nray), dtype=np.float32)
    cols = np.zeros((2 * Nside * Nproj * Nray), dtype=np.float32)
    vals = np.zeros((2 * Nside * Nproj * Nray), dtype=np.float32)
    idxend = 0

    for i in tqdm(range(0, Nproj)): # Loop over projection angles
        ang = angles[i] * np.pi / 180.
        # Points passed by rays at current angles
        xrayRotated = np.cos(ang) * offsets
        yrayRotated = np.sin(ang) * offsets
        xrayRotated[np.abs(xrayRotated) < 1e-8] = 0
        yrayRotated[np.abs(yrayRotated) < 1e-8] = 0

        a = -np.sin(ang)
        a = rmepsilon(a)
        b = np.cos(ang)
        b = rmepsilon(b)

        for j in range(0, Nray): # Loop rays in current projection
            #Ray: y = tx * x + intercept
            t_xgrid = (xgrid - xrayRotated[j]) / a
            y_xgrid = b * t_xgrid + yrayRotated[j]

            t_ygrid = (ygrid - yrayRotated[j]) / b
            x_ygrid = a * t_ygrid + xrayRotated[j]
            # Collect all points
            t_grid = np.append(t_xgrid, t_ygrid)
            xx = np.append(xgrid, x_ygrid)
            yy = np.append(y_xgrid, ygrid)
            # Sort the coordinates according to intersection time
            I = np.argsort(t_grid)
            xx = xx[I]
            yy = yy[I]

            # Get rid of points that are outside the image grid
            Ix = np.logical_and(xx >= -Nside / 2.0 * pixelWidth,
                                xx <= Nside / 2.0 * pixelWidth)
            Iy = np.logical_and(yy >= -Nside / 2.0 * pixelWidth,
                                yy <= Nside / 2.0 * pixelWidth)
            I = np.logical_and(Ix, Iy)
            xx = xx[I]
            yy = yy[I]

            # If the ray pass through the image grid
            if (xx.size != 0 and yy.size != 0):
                # Get rid of double counted points
                I = np.logical_and(np.abs(np.diff(xx)) <=
                                   1e-8, np.abs(np.diff(yy)) <= 1e-8)
                I2 = np.zeros(I.size + 1)
                I2[0:-1] = I
                xx = xx[np.logical_not(I2)]
                yy = yy[np.logical_not(I2)]

                # Calculate the length within the cell
                length = np.sqrt(np.diff(xx)**2 + np.diff(yy)**2)
                #Count number of cells the ray passes through
                numvals = length.size

                # Remove the rays that are on the boundary of the box in the
                # top or to the right of the image grid
                check1 = np.logical_and(b == 0, np.absolute(
                    yrayRotated[j] - Nside / 2 * pixelWidth) < 1e-15)
                check2 = np.logical_and(a == 0, np.absolute(
                    xrayRotated[j] - Nside / 2 * pixelWidth) < 1e-15)
                check = np.logical_not(np.logical_or(check1, check2))

                if np.logical_and(numvals > 0, check):
                    # Calculate corresponding indices in measurement matrix
                    # First, calculate the mid points coord. between two
                    # adjacent grid points
                    midpoints_x = rmepsilon(0.5 * (xx[0:-1] + xx[1:]))
                    midpoints_y = rmepsilon(0.5 * (yy[0:-1] + yy[1:]))
                    #Calculate the pixel index for mid points
                    pixelIndicex = ((np.floor(Nside / 2.0 - midpoints_y / pixelWidth)) * # noqa TODO reformat this
                        Nside + (np.floor(midpoints_x /
                        pixelWidth + Nside / 2.0)))
                    # Create the indices to store the values to the measurement
                    # matrix
                    idxstart = idxend
                    idxend = idxstart + numvals
                    idx = np.arange(idxstart, idxend)
                    # Store row numbers, column numbers and values
                    rows[idx] = i * Nray + j
                    cols[idx] = pixelIndicex
                    vals[idx] = length
            else:
                print(("Ray No.", j + 1, "at", angles[i],
                       "degree is out of image grid!"))
    # Truncate excess zeros.
    rows = rows[:idxend]
    cols = cols[:idxend]
    vals = vals[:idxend]
    A = np.array([rows, cols, vals], dtype=np.float32, order='C')
    return A


def rmepsilon(input):
    if (input.size > 1):
        input[np.abs(input) < 1e-10] = 0
    else:
        if np.abs(input) < 1e-10:
            input = 0
    return input

def load_data(vol_size, file_name):

    #sk-image loads tilt series as (z,y,x) so the axes need to be
    #swapped to return to (x,y,z)

    full_name = vol_size+file_name

    if full_name.endswith('.tiff'):
        tiltSeries = io.imread('Tilt_Series/'+ full_name)
        tiltSeries = np.array(tiltSeries, dtype=np.float32)
        tiltSeries = np.swapaxes(tiltSeries, 0, 2)
        file_name = file_name.replace('_tiltser.tiff', '')
    elif full_name.endswith('.tif'):
        tiltSeries = io.imread('Tilt_Series/'+ full_name)
        tiltSeries = np.array(tiltSeries, dtype=np.float32)
        tiltSeries = np.swapaxes(tiltSeries, 0, 2)
        file_name = file_name.replace('_tiltser.tif', '')
    elif full_name.endswith('.npy'):
        tiltSeries = np.load('Tilt_Series/' + full_name)
        file_name = file_name.replace('_tiltser.npy', '')

    return (file_name,tiltSeries)
    
