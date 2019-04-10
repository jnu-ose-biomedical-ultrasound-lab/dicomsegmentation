# mrReformat.py
# By John LaRocco
# Implemented for Jeju National University's Biomedical Ultrasound Lab (http://www.jejunu.ac.kr/eng/).

# 10 April 2019
# The purpose of this file is to run a Hessian filter over DICOM data to search for blood vessels.
# It searches each file, slice by slice.
# The slices are then integrated into a 3D matrix.
# It segments out blood vessel files based on HU values imported from the DICOM files, between 14 to 76.
# It finds all DICOM files in its home directory, and exports them as images and as a MAT file.

import cv2
from scipy import io
import pydicom as dicom
import os
import numpy
import scipy.misc
import numpy as np
from matplotlib import pyplot, cm
from skimage.filters import frangi, hessian
from skimage.feature import hessian_matrix, hessian_matrix_eigvals

PathDicom = "."
lstFiles = []
for dirName, subdirList, fileList in os.walk(PathDicom):
    for filename in fileList:

        if ".dcm" in filename.lower():
            lstFiles.append(os.path.join(dirName,filename))
def np2Mat(array,fn):
	data={'array':array}
	nameOut=str(fn)+'.mat'
	print(nameOut)
	io.savemat(nameOut,data)

def imageExport(array,fn):
    nameOut=str(fn)+'-nomask-veins.png'
    scipy.misc.imsave(nameOut, array)

def maskExport(array,fn):
    nameOut=str(fn)+'-mask.png'
    scipy.misc.imsave(nameOut, array)

def maskFinalExport(array,fn):
    nameOut=str(fn)+'-final-mask.png'
    scipy.misc.imsave(nameOut, array)

def aPrioriAppliedExport(array,fn):
    array[:,0:120]=0
    array[0:110]=0
    array[:,454::]=0
    array[475::]=0
    nameOut=str(fn)+'-apriori-veins.png'
    scipy.misc.imsave(nameOut, array)
    return(array)

def maskedVesselSegmentation(array):
    image, binary = boneSegmentation(array, [1, 300])
    kernel = np.ones((4,4),np.float32)/16
    img=np.squeeze(np.uint8(binary))
    img = cv2.filter2D(img,-1,kernel)
    img = cv2.erode(img, kernel, iterations=1)
    img = cv2.dilate(img, kernel, iterations=1)
    img = cv2.erode(img, kernel, iterations=1)
    img=frangi(img)
    hxx, hxy, hyy = hessian_matrix(img, sigma=3)
    i1, i2 = hessian_matrix_eigvals(hxx, hxy, hyy)
    img=i2/abs(i2.max())
    img=np.squeeze(np.uint8(img))
    threshold=img.max()-img.mean()
    img=np.squeeze(np.uint8(img))
    img = cv2.erode(img, kernel, iterations=1)
    img=np.squeeze(np.uint8(img))
    img = cv2.dilate(img, kernel, iterations=1)
    threshold=img.max()-img.mean()
    image = np.where((img > threshold), 255, 0)
    return(image)

def boneSegmentation(array, threshold):
    image = np.where((array > threshold[0]) & (array < threshold[1]), 255, 0)
    binary = np.where(image == 255, 1, 0)
    return image, binary

# Get ref file
RefDs = dicom.read_file(lstFiles[0])
# Load dimensions based on the number of rows, columns, and slices (along the Z axis)
ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFiles))
# Load spacing values (in mm)
ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))
x = numpy.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
y = numpy.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
z = numpy.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])
# The array is sized based on 'ConstPixelDims'
ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

# loop through all the DICOM files
for filenameDCM in lstFiles:

    dicomfile = dicom.read_file(filenameDCM)

    x=dicomfile.pixel_array
    xclone=np.copy(x)

    x=maskedVesselSegmentation(x)

    xclone[xclone >= 76] = 0
    xclone[xclone <= 14] = 0

    print(np.min(x))
    print(np.max(x))

    applyMask=xclone-x
    applyMask[applyMask < 0] = 0
    applyMask[applyMask > 75] = 75

    ArrayDicom[:, :, lstFiles.index(filenameDCM)] = x

    maskExport(applyMask,filenameDCM)


    imageExport(x,filenameDCM)
    final=aPrioriAppliedExport(x,filenameDCM)
    finalImage=final+applyMask
    maskFinalExport(finalImage,filenameDCM)
np2Mat(ArrayDicom,'rawmrifuseddicom')
