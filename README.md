For help, questions, corrections or feedback please contact: *awed2@cam.ac.uk*

# Analysing the effects of genotype on cortical network development using graph theory
This account is based on work that I completed during my Masters in the Neuronal Oscillations groups (Dept. of PDN, University of Cambridge) under the supervision of Dr Susanna Mierau and Prof Ole Paulsen.

*The final thesis will be uploaded on 18th September.*


## Scripts found in this repository were ultimately used for: 
(1) applying graph theory to cortical neurons cultured on microelectrode arrays; (2) comparing wild-type and Mecp2-deficient networks.


## Summary
These scripts take .mat files containing multielectrode array extracellular recordings an create spike matrices for each recording. That is, a sparse matrix is created where there is a binary vector for each channel indicating spike times (referred to as *spike trains*). Various analyses can be performed on spiking activity and one can also plot spikes and voltage traces in various ways. Bursting activity can then be detected using these spike matrices. This produces a cell variable which tells the user the burst start and end times as well as the electrodes containing spikes that contributed to the burst. Currently, functional connectivity is calculated by correlating spiking activity.

Currently, functional connectivity is calculated across entire spike trains rather than within bursts. This will be added in a future update.

## Steps

1. Record MEA activitiy using MC Rack software to produce .mcd files. Convert these to .raw files using the MC Data tool (available online). Then use MEAbatchConvert.m to convert the .raw files to .mat files. See the "mecp2" repository to get detailed instructions. The .mat files will contain a variable called "dat" — this contains a matrix with a row for each electrode and a column for each sample. Thus, if there are 60 electrodes in the array, there will be 60 row vectors containing voltage traces. There is also a vatiable called "channels" that tells the user which electrode each row corresponds to in terms of the electrode ID (e.g. electrode 78 is column eight, row 7 in the MEA). For example, if row 10 of *channels* is "78" then row 10 of *dat* is the voltage trace for the electrode in column 8, row 7 of the MEA.


2. Detect spikes using the .mat files. This creates spike matrices

3. Carry out spiking and bursting analyses on spike matrices

4. Create weighted adjacency matrices by correlating spike trains within spike matrices of each recording

5. Carry out functional connectivity and graph analyses.

**Steps 3–5 are done using network_features_MEA.m**
Currently, the script works on spike matrices with 60 channels. To exclude electrodes, set the spike train of that electrode to that of the reference electrode.

The network feature vector for each recording can be manually copied into excel. Desired feature can then be saved as a csv file and analysed with the R analysis script.

Future developments include further options for customisation; machine learning applied to the feature vector.
