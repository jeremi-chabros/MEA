# mphil
Scripts used for: (1) applying graph theory to cortical neurons cultured on microelectrode arrays; (2) comparing wild-type and Mecp2-deficient networks.

Use the network features script for the main analysis. First you need to get spike matrices. Currently, the script works on spike matrices with 60 channels. To exclude electrodes, set the spike train of that electrode to that of the reference electrode.

The network feature vector for each recording can be manually copied into excel. Desired feature can then be saved as a csv file and analysed with the R analysis script.

Future developments include further options for customisation; machine learning applied to the feature vector.
