# bloodHunt.py
# By John LaRocco
# Implemented for Jeju National University's Biomedical Ultrasound Lab (http://www.jejunu.ac.kr/eng/).

# 1 April 2019
# The purpose of this file is to run a Hessian filter over DICOM data to search for blood vessels.
# It searches each file, slice by slice.
# The slices are then integrated into a 3D matrix.
# A binary 3D matrix, based on the presence of possible blood vessels, is also output.
# TA similar search operation is performed for bone and hard tissue.
# Both of these operations result in a binary matrix. 
# In these matrices, 1 means that voxel has a blood vessel (or bone). 0 otherwise. 
# The first voxel coordinate for a suspected blood vessel will be automatically selected as a 'target.'
# Based on the 3D binary matrix, it can calculate a possible 'entry point' and angle, using a cost minimization algorithm.
# Currently, the body of the function is hard coded based on values of sample no100's DICOM files. 
# This function also exports the fused DICOM model and the blood and bone binary matrix as 3D .mat files (for MATLAB). 
# In each MAT file is a 3D matrix called 'array,' which contains the information from Python. 


from scipy import io
import os
import pydicom
import numpy as np
import pickle
from skimage.filters import frangi, hessian
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
from scipy.ndimage import rotate
from scipy import ndimage
from scipy import stats
from pylab import *
from scipy.ndimage import measurements,morphology
import matplotlib.patches as mpatches
from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb
from skimage.morphology import disk, dilation, binary_erosion, binary_closing
from skimage.filters import roberts, sobel
from scipy import ndimage as ndi
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D
from math import atan
import cv2




def load_patient(src_dir):
    slices = [pydicom.read_file(src_dir + '/' + s, force=True) for s in os.listdir(src_dir)]
#    slices.sort(key=lambda x: int(x.InstanceNumber))
    #slices.sort(key=lambda x: int(si for len(slices)))
    try:
        slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
    except:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)

    for s in slices:
        s.SliceThickness = slice_thickness
    return slices

def slicesInstanceNumber(slices):
    image = np.stack([s.pixel_array for s in slices])
    #image=frangi(image)
    return image

def boneSegmentation(array, threshold):

    image = np.where((array > threshold[0]) & (array < threshold[1]), 255, 0)
    binary = np.where(image == 255, 1, 0)
    return image, binary



def DICOMtoNumpyStack(slices):
    image=slicesInstanceNumber(slices)
    #image = np.stack([s.pixel_array for s in slices])
    image = image.astype(np.int16)
    image[image == -2000] = 0
    for slice_number in range(len(slices)):
        intercept = slices[slice_number].RescaleIntercept
        slope = slices[slice_number].RescaleSlope
        if slope != 1:
            image[slice_number] = slope * image[slice_number].astype(np.float64)
            image[slice_number] = image[slice_number].astype(np.int16)
        image[slice_number] += np.int16(intercept)

    return np.array(image, dtype=np.int16)




def DICOMtoNumpy(slices):
    image = np.stack([s.pixel_array for s in slices])
    #image=frangi(image)
    return image




def count_items(img1_f, it):
  im_open = morphology.binary_opening( \
    img1_f,ones((9,5)), iterations=it)
  labels_open, nbr_objects_open = measurements.label(im_open)
  return labels_open

def outputPNG(image, out_dir):
    for i in range(image.shape[0]):
        img_path = out_dir + "/img_" + str(i).rjust(4, '0') + ".png"
        cv2.imwrite(img_path, image[i])


def costMapGeneration(array):
# bone segmentation
    threshold = [1250, 3000]
    image1, binary1 = boneSegmentation(array, threshold)
# blood segmentation
    image, binary = boneSegmentation(array, [1, 250])
    img = np.uint8(image)

# vessel segmentation
    threshold2 = [1, 300]
    img=np.squeeze(img)
    kernel = np.ones((6,6),np.float32)/36
    img = cv2.filter2D(img,-1,kernel)
    imgB0=img
# best image
    img=np.squeeze(array)
    img = cv2.filter2D(img,-1,kernel)
    img0=img
    img = cv2.erode(img, kernel, iterations=2)
    img = cv2.dilate(img, kernel, iterations=2)
    img2 = cv2.erode(img, kernel, iterations=1)
    imf=img
    imf = cv2.dilate(imf, kernel, iterations=1)
    imf = cv2.erode(imf, kernel, iterations=1)
    imf = np.where((imf < 100), 0, 255)
    img11 = np.where((img < img.mean()), 0, 255)
    img21 = np.where((img2 < img.mean()), 0, 255)
    img3=img11+img21
    img4=np.where((img3 != 0), 255, 0)
    kernel = np.ones((2,2),np.float32)/4
    img = cv2.filter2D(img0,-1,kernel)
    hxx, hxy, hyy = hessian_matrix(img, sigma=1)
    i1, i2 = hessian_matrix_eigvals(hxx, hxy, hyy)
    img = cv2.erode(img, kernel, iterations=2)
    imgh = np.where((i2 < i2.mean()), 255, 0)
    img = cv2.dilate(img, kernel, iterations=2)
    img20 = cv2.erode(img, kernel, iterations=20)
    img20 = np.where((img20 < img20.mean()), 0, 255)
    otherPlace=np.where((img20 <= 0), 255, 0)
    otherPlace=np.where((otherPlace <= 0), 255, 0)
    img10 = cv2.erode(img, kernel, iterations=10)
    img10 = np.where((img10 < img10.mean()), 0, 255)
    otherPlace10=np.where((img10 <= 0), 255, 0)
    otherPlace10=np.where((otherPlace10 <= 0), 255, 0)
    img55 = cv2.erode(img, kernel, iterations=5)
    img55 = np.where((img55 < img55.mean()), 0, 255)
    otherPlace55=np.where((img55 <= 0), 255, 0)
    otherPlace55=np.where((otherPlace55 <= 0), 255, 0)
    img15 = cv2.erode(img, kernel, iterations=5)
    img15 = np.where((img15 < img15.mean()), 0, 255)
    otherPlace15=np.where((img15 <= 0), 255, 0)
    otherPlace15=np.where((otherPlace15 <= 0), 255, 0)
    OP=otherPlace15+otherPlace+otherPlace55+otherPlace10
    OP=np.where((OP > 0), 255, 0)
    img2 = cv2.erode(img, kernel, iterations=3)
    img = np.where((img < img.mean()), 0, 255)
    img2 = np.where((img2 < img.mean()), 0, 255)
    img5=img4-img2
    img6=np.where((img5 != 0), 255, 0)
    imgFrame=np.where((img20 <= img20.mean()), 0, 255)
    victorFrame=np.where((imf != 0), 0, 255)
    victorFrame=np.where((victorFrame <= 0), 0, 255)
    tangoZone=victorFrame+imgh

    bF=np.where((imgB0 < 255), 255, 0)
    OP1=OP-imf
    OP=np.where((OP <= 0), 255, 0)
    superZone=bF-OP
    superZone=np.where((superZone <= 0), 0, 255)
    superZone=superZone-img11
    superZone=np.where((superZone <= 0), 0, 255)
    superZone=superZone-OP1
    superZone=np.where((superZone <= 0), 0, 255)


    comboZone=tangoZone+img11-(binary)
    comboZone=np.where((comboZone <= 0), 255, 0)

    comboZone=np.where((comboZone <= 200), 255, 0)
    comboZone=np.where((comboZone <= 0), 255, 0)

    targetZone=comboZone
    comboZone=targetZone-otherPlace55
    comboZone=np.where((comboZone <= 0) & (comboZone < 255), 0, 255)
    binZone=superZone+comboZone
    binZone=np.where((binZone > 0), 255, 0)

    binV2=np.squeeze(binZone)
# bone cost
    binV1=np.squeeze(binary)
# all squares now have a cost
    costMapV=np.where(binV1==1,1000,1)
    costMapB=np.where(binV2==1,9000,1)

# final costmap
    costMap=np.squeeze(costMapV+costMapB)
    return(costMap,imgB0)


#set target equal to some nice negative number

def trajectoryPlannerSimple(costMap,target0,target1):
    A=costMap
    costMap[target0,target1]=-30000000
    topDown=sum(A[target0,0:target1]) #assuming directly above the target is 0 degrees of a circle
    bottomDown=sum(A[target0,target1:]) # assuming this is 180 degrees opposite of the 'top'
    leftEntrance=sum(A[0:target0,target1]) # 90 degrees, clockwise from top
    rightEntrance=sum(A[target0:,target1]) # 270 degrees, clockwise from top

    cost=np.array([topDown,leftEntrance,bottomDown,rightEntrance])
    angles=np.array([0,90,180,270])
    angleMin=angles[np.where(cost==cost.min())]
    return(angleMin,cost.min())


def trajectoryPlanner(costMap,target0,target1):
    A=costMap
    cost=np.zeros(359)
    cost[0]=sum(A[target0,0:target1]) #assuming directly above the target is 0 degrees of a circle
    for ii in range(0,359):
            costMap[target0,target1]=-30000000
            A=rotate(costMap,ii)
            cost[ii]=sum(A[target0,0:target1])
    angles=plt.imshow(np.squeeze(binZone), cmap = plt.get_cmap('gray'))

    angleMin=angles[np.where(cost==cost.min())]
    return(angleMin,cost.min())

def segmentAngle(costMap,target0,target1):

    angleMinSimple, costMinSimple=trajectoryPlannerSimple(costMap,target0,target1)

    angleMin, costMin=trajectoryPlanner(costMap,target0,target1)

    recAngle=np.array([angleMinSimple,angleMin])
    costCompare=np.array([costMinSimple,costMin])
    finalCostMin=costCompare.min()
    finalAngle=recAngle[np.where(costCompare==finalCostMin)]
    return(finalAngle,finalC5ostMin)


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


def trajectoryPlanner3D(costMap,target0,target1,target2):
    A=costMap
    A[target0,target1,target2]=-30000000
    distance=360*3
    xx=np.array(np.zeros(360))
    yy=xx
    zz=yy
    costX0=sum(A[0:target0,target1,target2])
    costY0=sum(A[target0,0:target1,target2])
    costZ0=sum(A[target0,target1,0:target2])
    for xi in range(0,359):
        xplane=np.squeeze(A[target0,:,:])
        x2=rotate(xplane,xi)
        xx[xi]=sum(x2[0,0:target1])+sum(x2[1,0:target2])
    costX=xx.min()
    xAng=np.where((xx==xx.min()) & (xx!=0))
    for yi in range(0,359):
        yplane=np.squeeze(A[:,target1,:])
        y2=rotate(yplane,yi)
        yy[yi]=sum(y2[0,0:target0])+sum(y2[1,0:target2])
    costY=yy.min()
    yAng=np.where((yy==yy.min()) & (yy!=0))

    for zi in range(0,359):
        zplane=np.squeeze(A[:,:,target2])
        z2=rotate(zplane,zi)
        zz[zi]=sum(z2[0,0:target0])+sum(z2[1,0:target1])
    costZ=zz.min()
    zAng=np.where((zz==zz.min()) & (zz!=0))
    minCost=np.array([costX,costY,costZ])
    if costX > costX0:
        xAng=0
        minCost[0]=costX0
    if costY > costY0:
        yAng=0
        minCost[1]=costY0
    if costZ > costZ0:
        zAng=0
        minCost[2]=costZ0n
    totalCost=sum(minCost)
    xPro=minCost[0]/totalCost
    yPro=minCost[1]/totalCost
    zPro=minCost[2]/totalCost
    outVec=np.zeros(3)
    outVec[0]=np.asarray(xAng**xPro)
    outVec[1]=np.asarray(yAng**yPro)
    outVec[2]=np.asarray(zAng**zPro)
    return(xAng,yAng,zAng,minCost,outVec)


def wheelCost(costMap,target0,target1):
    # Calculate cost for each 45 degree increment.
    a1=np.trace(costMap[:target0,:target1])  # 1 + 5: 315

    a2=np.sum(costMap[:target0,:target1], axis=0)[1]  # 2 + 5: 0

    a3=np.trace(np.rot90(costMap[:target0,(target1-1):]))  # 5 + 3: 45

    a4=np.sum(costMap[:target0,(target1-1):], axis=1)[1]  # 5 + 6: 90

    a5=np.trace(costMap[(target0-1):,(target1-1):])  # 5 + 9: 125

    a6=np.sum(costMap[(target0-1):,(target1-1):], axis=0)[0]  # 5 + 8: 180

    a7=np.trace(np.rot90(costMap[(target0-1):,:target1]))  # 7 + 5: 225

    a8=np.sum(costMap[(target0-1):,:target1], axis=1)[0]  # 4 + 5: 270
    a=np.array([a1,a2,a3,a4,a5,a6,a7,a8])
    anglePos=np.where(a==a.min())
    return(anglePos,a)


def maskFusion(img,mask):
    # Combine the masks.
    fused=img-mask

    img=np.where((fused >= 0), 255, 0)
    outImg=np.where((img == 0), 255, 0)
    return(outImg)



def costSlicer(outImg,binary):
    # fuse the cost maps of blood and bone
    if abs(outImg.max()==0):
        outImg=outImg
    else:
        outImg=outImg/abs(outImg.max())

    bloodImg=np.uint8(outImg)

    bloodImg=np.where(bloodImg>=1,9000,0)

    if abs(binary.max()==0):
        boneImg=np.uint8(abs(binary/binary.max()))
    else:
        boneImg=np.uint8(abs(binary))
    boneImg=np.where(boneImg>=1,6000,0)
    costSlice=np.squeeze(bloodImg+boneImg)
    return(costSlice)


def starPoint(anglePos1):
    degree=anglePos1[0].max()
    if (degree >= 0) & (degree <= 30):
        cd = 0
    if (degree >> 30) & (degree <= 60):
        cd = 1
    if (degree >> 60) & (degree <= 105):
        cd = 2
    if (degree >> 105) & (degree <= 150):
        cd = 3
    if (degree >> 150) & (degree <= 195):
        cd = 4
    if (degree >> 195) & (degree <= 240):
        cd = 5
    if (degree >> 240) & (degree <= 285):
        cd = 6
    if (degree >> 285) & (degree <= 330):
        cd = 7
    if (degree >> 330) & (degree <= 359):
        cd = 0
    if cd==0: # val
        dx=0
        dy=-1

    if cd==1: # val
        dx=1
        dy=-1

    if cd==2:# val
        dx=1
        dy=0

    if cd==3:# val
        dx=1
        dy=1

    if cd==4: #val
        dx=0
        dy=1

    if cd==5:# val
        dx=-1
        dy=1

    if cd==6:# val
        dx=-1
        dy=0

    if cd==7:# val
        dx=-1
        dy=-1
    return(dx,dy,cd)


def endPointCalculator(anglePos,dimsize1,dimsize2,targetA,targetB):
    cord1=0
    cord2=0
# Note: This values require debugging and testing.
    if anglePos==0:
        cord1=0
        cord2=(dimsize2-1)

    if anglePos==1:
        cord1=targetA
        cord2=(dimsize2-1)
    if anglePos==2:
        cord1=(dimsize1-1)
        cord2=(dimsize2-1)
    if anglePos==3:
        cord1=(dimsize1-1)
        cord2=targetB
    if anglePos==4:
        cord1=(dimsize1-1)
        cord2=0
    if anglePos==5:
        cord1=targetA
        cord2=0
    if anglePos==6:
        cord1=0
        cord2=0
    if anglePos==7:
        cord1=0
        cord2=targetB

    return(cord1,cord2)


def sphereCost(costMap,target0,target1,target2):
# This function was intended to calculate 26 possible trajectors for the input point.
# First, you size length of each dimension in the 3D array.
# Then, wheelCost is used to calculate the lowest cost angle (anglePos) and cost value vector (a) for each plane.
# Following this, the lowest cost trajectory is planned. The end point of each trajectory is calculated.
# All values are compared to find the lowest total cost trajectories.
# The cost map values in between the target and end point values will be converted to a high contrast value to allow for easier plotting.


    costMap=np.asarray(costMap)

    xsize=len(np.squeeze(costMap[:,0,0]))

    ysize=len(np.squeeze(costMap[0,:,0]))

    zsize=len(np.squeeze(costMap[0,0,:]))
    costMap=costMap+1
    direct1=np.sum(costMap[0:target0,0:target1,0:target2])

    direct2=np.sum(costMap[target0:(xsize-1),target1:(ysize-1),target2:(zsize-1)])

    axisx1=np.sum(costMap[0:target0,0,0])
    axisx1=np.sum(costMap[0:target0,target1,target2])
    axisx2=np.sum(costMap[target0:(xsize-1),0,0])
    axisx2=np.sum(costMap[target0:(xsize-1),target1,target2])
    minX=min([axisx1,axisx2])

    axisy1=np.sum(costMap[0,0:target1,0])
    axisy1=np.sum(costMap[target0,0:target1,target2])
    axisy2=np.sum(costMap[0,target1:(ysize-1),0])
    axisy2=np.sum(costMap[target0,target1:(ysize-1),target2])
    minY=min([axisy1,axisy2])

    axisz1=np.sum(costMap[0,0,0:target2])
    axisz1=np.sum(costMap[target0,target1,0:target2])
    axisz2=np.sum(costMap[0,0,target2:(zsize-1)])
    axisz2=np.sum(costMap[target0,target1,target2:(zsize-1)])
    minZ=min([axisz1,axisz2])

    defFrame=min([axisx1,axisx2,axisy1,axisy2,axisz1,axisz2,direct1,direct2])

    if (defFrame==axisz1) or (defFrame==axisy1) or (defFrame==axisz1) or (defFrame==direct1):
        dx=0
        dy=0
        dz=0

    if defFrame==axisz2:
        dx=0
        dy=0
        dz=(zsize-1)

    if defFrame==axisy2:
        dx=0
        dy=(ysize-1)
        dz=0

    if defFrame==axisx2:
        dx=(xsize-1)
        dy=0
        dz=0

    if defFrame==direct2:
        dx=(xsize-1)
        dy=(ysize-1)
        dz=(zsize-1)




    # Here you would use the wheelCost algorithm for each plane/
    #XY
    anglePos1,a1=wheelCost(np.squeeze(costMap[:,:,target2]),target0,target1)
    tracA=np.sum(costMap[0:target0,0:target1,target2])
    tracB=np.sum(costMap[0:target0,target1,0:target2])
    tracC=np.sum(costMap[target0,0:target1,0:target2])

    #XZ
    anglePos2,a2=wheelCost(np.squeeze(costMap[:,target1,:]),target0,target2)
    anglePos3,a3=wheelCost(np.squeeze(costMap[target0,:,:]),target1,target2)

    cordA1,cordA2=endPointCalculator(anglePos1,xsize,ysize,target0,target1)
    x=np.sort([target0,cordA1])
    y=np.sort([target1,cordA2])
    trajA=np.sum(np.squeeze(costMap[x[0]:x[1],y[0]:y[1],target2]))


    cordB1,cordB2=endPointCalculator(anglePos2,xsize,zsize,target0,target2)
    x=np.sort([target0,cordB1])
    z=np.sort([target2,cordB2])
    trajB=np.sum(np.squeeze(costMap[x[0]:x[1],target1,z[0]:z[1]]))
    trajAB=np.sum(np.squeeze(costMap[x[0]:x[1],y[0]:y[1],target2]))
    cordC1,cordC2=endPointCalculator(anglePos3,ysize,zsize,target1,target2)
    y=np.sort([target1,cordC1])
    z=np.sort([target2,cordC2])
    trajC=np.sum(np.squeeze(costMap[target0,y[0]:y[1],z[0]:z[1]]))

    defTraj=min([trajA, trajAB, trajB, trajC])


    entryCoord=np.array([cordA1,cordA2,cordB2])
    if defTraj == trajA:
        dx1=cordA1
        dy1=cordA2
        dz1=target2
    if defTraj == trajAB:
        dx1=cordB1
        dy1=cordA2
        dz1=target2
    if defTraj == trajB:
        dx1=cordB1
        dy1=target1
        dz1=cordB2
    if defTraj == trajC:
        dx1=target0
        dy1=cordC1
        dz1=cordC2
    if (abs(defTraj) < abs(defFrame)):
        defFrame=defTraj
        dx=dx1
        dy=dy1
        dz=dz1


    entryCoord=np.array([dx,dy,dz])
    x=np.sort([target0,entryCoord[0]])
    y=np.sort([target1,entryCoord[1]])
    z=np.sort([target2,entryCoord[2]])


    maxVal=100*costMap.max()

    costAp=np.copy(costMap)

    costAp[x[0]:x[1],y[0]:y[1],z[0]:z[1]]=maxVal
    costAp=costAp/maxVal

    return(costAp,x,y,z,entryCoord,maxVal,anglePos1,anglePos2)



def rawDicomMiner(file_name,filename,range1):
# Function exports the DICOM values without any segmentation. 
    nSlices=len(range1)-1
    slicer = [pydicom.read_file(filename, force=True)]
    image = DICOMtoNumpy(slicer)

    imgsize=image.shape
    imgsize=np.asarray(imgsize)
    xBor=imgsize[1]
    yBor=imgsize[2]
    costMap=np.zeros([xBor,yBor,nSlices])

    for i in range (0, nSlices):
        filename=file_name+str(range1[i])+'.dcm'
      #  filename=file_name+str(range1[i])
        slices = [pydicom.read_file(filename, force=True)]
        image = DICOMtoNumpy(slices)

        if i==0:
            ind=0
            newOrder=ind
            costMap[:,:,newOrder]=np.squeeze(image)

        if i >>1:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)

        if i==nSlices:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)
	return(costMap)


def rawDicomLister(resultList):
# Function exports the DICOM values without any segmentation. 
    nSlices=len(resultList)-1
    slicer = [pydicom.read_file(resultList[0], force=True)]
    image = DICOMtoNumpy(slicer)

    imgsize=image.shape
    imgsize=np.asarray(imgsize)
    xBor=imgsize[1]
    yBor=imgsize[2]
    costMap=np.zeros([xBor,yBor,nSlices])

    for i in range (0, nSlices):
        #filename=file_name+str(range1[i])+'.dcm'
      #  filename=file_name+str(range1[i])
        filename=resultList[i]
        slices = [pydicom.read_file(filename, force=True)]
        image = DICOMtoNumpy(slices)

        if i==0:
            ind=0
            newOrder=ind
            costMap[:,:,newOrder]=np.squeeze(image)

        if i >>1:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)

        if i==nSlices:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)
	return(costMap)




def superDicomMiner(resultList):
# Function exports the DICOM values with any segmentation into tissue types: regular, bone, and soft (blood).
    nSlices=len(resultList)-1
    slicer = [pydicom.read_file(resultList[0], force=True)]
    image = DICOMtoNumpy(slicer)

    imgsize=image.shape
    imgsize=np.asarray(imgsize)
    xBor=imgsize[1]
    yBor=imgsize[2]
    costMap=np.zeros([xBor,yBor,nSlices])
    boneMap=np.zeros([xBor,yBor,nSlices])
    bloodMap=np.zeros([xBor,yBor,nSlices])
	
    for i in range (0, nSlices):
        filename=resultList[i]
        slices = [pydicom.read_file(filename, force=True)]
        image = DICOMtoNumpy(slices)
	
    
        imageBone, binaryBone = boneSegmentation(image, [300, 3000])
        imageBone = np.uint8(binaryBone)
        imageBlood, binaryBlood = boneSegmentation(image, [1, 250])
        imageBlood = np.uint8(binaryBlood)

        if i==0:
            ind=0
            newOrder=ind
            costMap[:,:,newOrder]=np.squeeze(image)
            boneMap[:,:,newOrder]=np.squeeze(imageBone)
            bloodMap[:,:,newOrder]=np.squeeze(imageBlood)

        if i >>1:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)
            boneMap[:,:,newOrder]=np.squeeze(imageBone)
            bloodMap[:,:,newOrder]=np.squeeze(imageBlood)
			
        if i==nSlices:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)
            boneMap[:,:,newOrder]=np.squeeze(imageBone)
            bloodMap[:,:,newOrder]=np.squeeze(imageBlood)



	return(costMap,boneMap,bloodMap)







def dicomMiner(file_name,filename,range1):
# Function exports the DICOM values with any segmentation into tissue types: regular, bone, and soft (blood).
    nSlices=len(range1)-1
    slicer = [pydicom.read_file(filename, force=True)]
    image = DICOMtoNumpy(slicer)

    imgsize=image.shape
    imgsize=np.asarray(imgsize)
    xBor=imgsize[1]
    yBor=imgsize[2]
    costMap=np.zeros([xBor,yBor,nSlices])
    boneMap=np.zeros([xBor,yBor,nSlices])
    bloodMap=np.zeros([xBor,yBor,nSlices])
	
    for i in range (0, nSlices):
        filename=file_name+str(range1[i])+'.dcm'
        slices = [pydicom.read_file(filename, force=True)]
        image = DICOMtoNumpy(slices)
        
        imageBone, binaryBone = boneSegmentation(image, [300, 3000])
        imageBone = np.uint8(binaryBone)
        imageBlood, binaryBlood = boneSegmentation(image, [1, 250])
        imageBlood = np.uint8(binaryBlood)

        if i==0:
            ind=0
            newOrder=ind
            costMap[:,:,newOrder]=np.squeeze(image)
            boneMap[:,:,newOrder]=np.squeeze(imageBone)
            bloodMap[:,:,newOrder]=np.squeeze(imageBlood)

        if i >>1:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)
            boneMap[:,:,newOrder]=np.squeeze(imageBone)
            bloodMap[:,:,newOrder]=np.squeeze(imageBlood)
			
        if i==nSlices:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(image)
            boneMap[:,:,newOrder]=np.squeeze(imageBone)
            bloodMap[:,:,newOrder]=np.squeeze(imageBlood)



	return(costMap,boneMap,bloodMap)


def np2Mat(array,fn):
	data={'array':array}
	nameOut=str(fn)+'.mat'
	print(nameOut)
	io.savemat(nameOut,data)


def veinMiner(file_name,filename,range1):
    nSlices=len(range1)-1
    slicer = [pydicom.read_file(filename, force=True)]
    image = DICOMtoNumpy(slicer)

    imgsize=image.shape
    imgsize=np.asarray(imgsize)
    xBor=imgsize[1]
    yBor=imgsize[2]
    costMap=np.zeros([xBor,yBor,nSlices])

    for i in range (0, nSlices):

        filename=file_name+str(range1[i])+'.dcm'
        slices = [pydicom.read_file(filename, force=True)]
        image = DICOMtoNumpy(slices)
        filenameOut1='unaltered'+str(file_name)+str(range1[i])+'.png'
        imageV=np.squeeze(np.asarray(image))

        imsave(filenameOut1,imageV)
        imgsize=image.shape
        imgsize=np.asarray(imgsize)
        xBor=imgsize[1]
        yBor=imgsize[2]

#plt.imshow(image, cmap = plt.get_cmap('gray'))
#plt.show()
        imageVeins=maskedVesselSegmentation(image)
        imageVeins=np.where((imageVeins == 0), 255, 0)
        imageVeins=imageVeins/imageVeins.max()
        imageVeins=np.where((imageVeins == 1), 0, 1)

        if i==0:
            ind=0
            newOrder=ind
            costMap[:,:,newOrder]=np.squeeze(imageVeins)

        if i >>1:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(imageVeins)

        if i==nSlices:
            ind=i-1
            newOrder=int(ind)
            costMap[:,:,newOrder]=np.squeeze(imageVeins)
        filenameOut=str(file_name)+str(range1[i])+'.png'
        imsave(filenameOut,imageVeins)
        #filenameOut1=str('unaltered')+str(file_name)+str(range1[i])+'.png'
        #imsave(filenameOut1,image)
    return(costMap)

# Main test
target0=1
if target0==0:
    target0=1
target1=11
if target1==0:
    target1=1
target2=3
if target2==0:
    target2=1
targets=np.array([target0,target1,target2])

lbnd=0
ubnd=284
results=[]

for f in os.listdir('.'):
	if f.endswith('.dcm'):
		results.append(f)



#file_name = 'image'
#range1=range(lbnd,ubnd)
#filename=file_name+str(range1[0])+'.dcm'
#filename=file_name+str(range1[0])
#rawDicomMap=rawDicomMiner(file_name,filename,range1)
rawDicomMap=rawDicomLister(results)
np2Mat(rawDicomMap,'rawdicom')
dicomMap,boneMap,bloodMap=superDicomMiner(results)
#dicomMap,boneMap,bloodMap=dicomMiner(file_name,filename,range1)
np2Mat(dicomMap,'dicom')
np2Mat(boneMap,'bone')
np2Mat(bloodMap,'blood')
