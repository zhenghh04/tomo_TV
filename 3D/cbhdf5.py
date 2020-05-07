#!/usr/bin/env python
import numpy as np
import h5py
import argparse
from Utils.pytvlib import load_data
parser = argparse.ArgumentParser(description='Combine data')
parser.add_argument("--tiltAngles", '-a', default='au_sto_tiltAngles.npy')
parser.add_argument("--tiltSeries", '-s', default='1024_au_sto_tiltser.npy')
parser.add_argument("--output", '-o', default='output.h5')
args = parser.parse_args()
series=np.load(args.tiltSeries)
print(series)
h5 = h5py.File(args.output, 'w')
h5.create_dataset('tiltAngles', data=np.load(args.tiltAngles))
fs=h5.create_dataset('tiltSeries', data=series)
h5.close()

