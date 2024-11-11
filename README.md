# Dependence of microstructure classification accuracy on crystallographic data representation
GitHub respository for those interested in training models with specific crystallographic data representations. 

![alt text](https://github.com/Shrunalp/cnn-mat/blob/main/sample_splicing.jpg?raw=true#center)

## Abstract
Convolutional neural networks are increasingly being used to analyze and classify material microstructures, motivated
by the possibility that they will be able to identify relevant microstructural features more efficiently and impartially
than human experts. While up to now convolutional neural networks have mostly been applied to light optimal mi-
croscopy and scanning electron microscope micrographs, application to EBSD micrographs will be increasingly common
as rational design generates materials with unknown textures and phase compositions. This raises the question of how
crystallographic orientation should be represented in such a convolutional neural network, and whether this choice has
a significant effect on the network’s analysis and classification accuracy. Four representations of orientation information
are examined and are used with convolutional neural networks to classify five synthetic microstructures with varying
textures and grain geometries. Of these, a spectral embedding of crystallographic orientations in a space that respects
the crystallographic symmetries performs by far the best, even when the network is trained on small volumes of data
such as could be accessible by practical experiments.

## Installation

Clone the Github repoistory:
```
git clone https://github.com/Shrunalp/cnn-mat.git
```
or using Git CLI 
```
gh repo clone Shrunalp/cnn-mat
```


## Converting DREAM.3D data into varying GVD block representations 
This implementation requires that the FID/IPF/PQ data generated from DREAM.3D is a .txt file extension. 
Use ```txt_to_np``` to convert the txt file to a numpy array. The function ```np_to_block``` will then convert 
the corresponding numpy array into a GVD block representation. 

Converting the PQ data into the SEQ data requires using the MATLAB. Simply go the MATLAB file ```test_pq_embedding```
and change the file name to the file name of the PQ .txt data from DREAM.3D. The implementation from our paper uses
embedding ```n_dim =4``` but this can be varied according to your preference. The data will processed and returned as 
a SEQ txt file which can be used with the instructions above to create the block representation.

## Hyperparameter sweep for training CNNs
(In progress)

## Authors
**PIs**
* Dr. Benjamin Schweinhart
* Dr. Jerry Mason
* Dr. Tyrus Berry

**Graduate Students**
* Dylan Miley
* Shrunal Pothagoni
