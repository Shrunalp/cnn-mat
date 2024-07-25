import csv
import numpy as np
import random
import tensorflow as tf

def txt_to_np(path, IPF_files, Quat_files):

# Converts DREAM3D txt files to compatible numpy files (.npy)
# 
# Inputs:
#   path        - Specify the location of the txt files in your directory as a string.
# 
#   IPF_files    - Names of the IPF/FID files as a list of strings.
#
#   Quat_files   - Names of the Grid data files as a list of strings.
# 
# Outputs:
#   Numpy Array - Returns a reformated numpy array for both the FID/IPF and PQ datasets.
#


### Open Data ###

    FID_IPF = []
    Phase_Quat = []

    for i in range(len(IPF_files)):
        data = open(path+IPF_files[i], "r")
        FID_IPF.append(data)

        data = open(path+Quat_files[i], "r")
        Phase_Quat.append(data)

### Convert FID/IPF to Compatible Data ###

    FID_IPF_Data = []
    matl = []

    for i in range(len(IPF_files)):
        for j in FID_IPF[i]:
            matl.append(j.replace("\n","").split(','))
        FID_IPF_Data.append(matl)
        matl = []

    FID_IPF_Data = np.array(FID_IPF_Data, dtype=float)
    name_FID = input("What would you like to name the FID/IPD data?")
    np.save(name_FID, FID_IPF_Data)

### Convert PQ to Compatible Data ###

    PQ_Data = []
    for i in range(len(Quat_files)):
        matl = []
        for j in Phase_Quat[i]:
            matl.append(j.replace("\n","").split(','))
        PQ_Data.append(matl)

    PQ_Data = np.array(PQ_Data, dtype=float)
    name_PQ = input("What would you like to name the PQ data?")
    np.save(name_PQ, PQ_Data)


#####################################################################################################


def np_to_block(FID_IPF_Data, Phase_Quat_Data, Block_Size=(4,512,512,512)):

# Converts DREAM3D numpy files (.npy) to block data for 
# 
# Inputs:
#   FID_IPF_Data        - Name of the IPF/FID file as a string.
# 
#   Phase_Quat_Data     - Names of the PQ file as a string.
#
#   Block_Size          - Shape of grid data as a tuple of ints.
# 
# Outputs:
#   Numpy Array - Returns a reformated numpy array block given by Block_Size
#


### Reformat Data for Training ###

### FID 2D Array ###
    img_data_FID = np.zeros(Block_Size)
    for i in range(Block_Size[1]):
        for j in range(Block_Size[2]):
            for k in range(Block_Size[3]):
                img_data_FID[0][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][0]
                img_data_FID[1][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][0]
                img_data_FID[2][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][0]
                img_data_FID[3][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][0]

### RGB 2D Array using IPF Colors ###
    img_data_IPF = np.zeros(Block_Size)
    for i in range(Block_Size[1]):
        for j in range(Block_Size[2]):
            for k in range(Block_Size[3]):
                img_data_IPF[0][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][1]
                img_data_IPF[1][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][2]
                img_data_IPF[2][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][3]
                img_data_IPF[3][i][j][k] = FID_IPF_Data[i+Block_Size[2]*j+(Block_Size[1]*Block_Size[2])*k + 1][1]

### CMYB 2D Array using Quaternions ###
    img_data_PQ = np.zeros(Block_Size)
    for i in range(Block_Size[1]):
        for j in range(Block_Size[2]):
            for k in range(Block_Size[3]):
                x = int(img_data_FID[0][i][j][k])
                img_data_PQ[0][i][j][k] = Phase_Quat_Data[x][1]
                img_data_PQ[1][i][j][k] = Phase_Quat_Data[x][2]
                img_data_PQ[2][i][j][k] = Phase_Quat_Data[x][3]
                img_data_PQ[3][i][j][k] = Phase_Quat_Data[x][4]

    print('Preprocessing is finished.')

### Naming ###
    name_FID = input("What would you like to name the FID Block?")
    np.save(name_FID, img_data_FID)

    name_IPF = input("What would you like to name the IPF Block?")
    np.save(name_IPF, img_data_IPF)

    name_PQ = input("What would you like to name the PQ Block?")
    np.save(name_PQ , img_data_PQ)


#####################################################################################################


def b2w(Block_data):

# Converts DREAM3D numpy block data (.npy) to imgs and label tensors for Tensorflow
# 
# Inputs:
#   Block_data  - Name of the block data as strings.
# 
# Outputs:
#   List - Returns a list with imgs and label tensors
#
    matl = [np.load(i).T for i in Block_data]

### Split Blocks into Windows of size 64x64 over a 10p increment ###
    def convert_to_window(block_matl):
        wind = []
        for i in range(51):
            for j in range(8):
                for k in range(8):
                    x = block_matl[i*10, 64*j:64*(j+1), 64*k:64*(k+1), 0:4]
                    wind.append(x)
        wind = np.array(wind)

        return wind

    window_data = np.array([convert_to_window(i) for i in matl])
    lb = np.array([np.full((len(window_data[0]),1), int(i)) for i in range(len(matl))])

    imgs = np.vstack(tuple(window_data))
    labels = np.vstack(tuple(lb))

### Combine the two lists using zip ###
    combined = list(zip(imgs, labels))

### Shuffle the combined list ###
    random.shuffle(combined)

### Unpack the shuffled list into separate lists ###
    imgs, labels = zip(*combined)

### Convert to tensors ###
    imgs = tf.convert_to_tensor(np.array(imgs))
    labels = tf.one_hot(np.array(labels), depth=len(matl))
    labels =  tf.squeeze(labels, axis=1)

    return [imgs,labels]



