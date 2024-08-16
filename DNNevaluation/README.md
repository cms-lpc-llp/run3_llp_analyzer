# Add DNN branches to the analysis ntuples

### Train and store the .h5 model

April 17th file: training_CA0p6_NoMerging_FlatClusterSize0to500_230415.h5

- Input ntuples: CA clustering with dR = 0.6, no merging. 

- Trained with a ClusterSize flat distribution from 50 to 500. Balanced classes (10k signal DT clusters, 10k signal CSC clusters, 10k bkg DT clusters, 10k bkg CSC clusters).

- Training variables: Cluster_isCSC, Cluster_XSpread, Cluster_YSpread, Cluster_ZSpread, Cluster_XYSpread, Cluster_RSpread, frac_s1, frac_s2, frac_s3, frac_s4, frac_rw1, frac_rw2, frac_rw3

- DNN classifier architecture: 

fully-connected with hidden layers (64, 64, 32, 16, 8)
binary cross-entropy as loss function
Adam as optimizer 
learning rate = 0.0005
number of epochs = 1000 
10% of the training dataset used as validation dataset
Early stopping with pacient = 30 epochs on the loss of the validation dataset to control the overtraining. Early stopped on epoch 216 -> loss: 0.2271 - accuracy: 0.9076 - val_loss: 0.3398 - val_accuracy: 0.8510

May 2024 file: training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5

- Adversarial DNN trained with MC + beam halo data clusters
- More details: https://indico.cern.ch/event/1391693/#2-machine-learning-studies-for

### Get list of files

python3 make_list_of_files.py

### submit jobs

condor_submit submit.sub


