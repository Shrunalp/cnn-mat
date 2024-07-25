import csv
import numpy as np
import random

### combines fid_ipf and quaternion datasets. takes filepath and filename as input
# saves ndarrays of floats for assembled fid_ipf and quaternion datasets.
# IMPORTANT to have coreresponding ipf and quat data
# in corresponding elements of IPF_files and Quat_files 

def concat_and_np(path, IPF_files, Quat_files):

    FID_IPF = []
    Phase_Quat = []

    for i in range(len(IPF_files)):
        data = open(path+IPF_files[i], "r")
        FID_IPF.append(data)

        data = open(path+Quat_files[i], "r")
        Phase_Quat.append(data)

### FID_IPF ###

    matl = []
    FID_IPF_Data = []
    for i in range(len(IPF_files)):
        for j in FID_IPF[i]:
            matl.append(j.replace("\n","").split(','))
        FID_IPF_Data.append(matl)
        matl = []

### deletes unnessary column labels and converts remaining to ndarray floats ### 
    fid1 = np.delete(FID_IPF_Data[0],0,0)
    fid2 = np.delete(FID_IPF_Data[1],0,0)
    fid3 = np.delete(FID_IPF_Data[2],0,0)
    fid4 = np.delete(FID_IPF_Data[3],0,0)

    fid1 = np.array(fid1, dtype=float)
    fid2 = np.array(fid2, dtype=float)
    fid3 = np.array(fid3, dtype=float)
    fid4 = np.array(fid4, dtype=float)

### finds max value (inherently the highest FID and thus the number of grains ###
    max1 = np.amax(fid1)
    max2 = np.amax(fid2,axis=0)[0]
    max3 = np.amax(fid3,axis=0)[0]
    max4 = np.amax(fid4,axis=0)[0]

### adds max fid of X1 to fids of X2, adds max of X1 and X2 to X3 etc. ###
    for i in range(fid2.shape[0]):
        fid2[i][0] += max1
    for i in range(fid3.shape[0]):
        fid3[i][0] += max1+max2
    for i in range(fid4.shape[0]):
        fid4[i][0] += max1+max2+max3

### stitches everything together ###
    fid_cat = np.concatenate((fid1,fid2,fid3,fid4), axis=0)
    fid_cat = np.array(fid_cat,dtype=float)
    print('concatenated fid:',fid_cat.shape)
    print('max fid val:',np.amax(fid_cat))

    name_FID = input("What would you like to name the FID/IPD data?")
    np.save(name_FID, fid_cat)

### PQ ###

    PQ_Data = []
    for i in range(len(Quat_files)):
        matl = []
        for j in Phase_Quat[i]:
            matl.append(j.replace("\n","").split(','))
        PQ_Data.append(matl)

    pq1 = (np.delete(PQ_Data[0],[0,1],0))
    pq2 = (np.delete(PQ_Data[1],[0,1],0))
    pq3 = (np.delete(PQ_Data[2],[0,1],0))
    pq4 = (np.delete(PQ_Data[3],[0,1],0))

    pq1 = np.array(pq1,dtype=float)
    pq2 = np.array(pq2,dtype=float)
    pq3 = np.array(pq3,dtype=float)
    pq4 = np.array(pq4,dtype=float)

    pq_cat = np.vstack((pq1,pq2,pq3,pq4))
    print('new shape for PQ:',pq_cat.shape)

    name_PQ = input("What would you like to name the PQ data?")
    np.save(name_PQ, pq_cat)

############################################################################################
