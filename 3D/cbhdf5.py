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
print(np.load(args.tiltSeries, allow_pickle=True))
series=np.load(args.tiltSeries, allow_pickle=True)[1]
h5 = h5py.File(args.output, 'w')
h5.create_dataset('tiltAngles', data=np.load(args.tiltAngles, allow_pickle=True))
fs=h5.create_dataset('tiltSeries', data=series)
h5.close()

