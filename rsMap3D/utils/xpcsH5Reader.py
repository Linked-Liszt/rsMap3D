'''
 Copyright (c) 2022, UChicago Argonne, LLC
 See LICENSE file.
'''
'''
Methods to assist with the reading of various H5 metadata files. 
'''
import h5py
import numpy as np

def getDetectorROI(metaFile):
    with h5py.File(metaFile, 'r') as dataF:
        return [
                dataF['measurement/instrument/detector/roi/y1'][0,0],
                dataF['measurement/instrument/detector/roi/y2'][0,0] + 1, 
                dataF['measurement/instrument/detector/roi/x1'][0,0],
                dataF['measurement/instrument/detector/roi/x2'][0,0] + 1
            ]


def getImages(dataFile, startImage, numImages):
    with h5py.File(dataFile, 'r') as dataF:
        return np.asarray(dataF['entry/instrument/detector/data'][startImage:startImage+numImages,:,:])


def getNumImages(dataFile):
    with h5py.File(dataFile, 'r') as dataF:
        return dataF['entry/instrument/detector/data'].shape[0]