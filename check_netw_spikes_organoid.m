% check network spikes in organoids

% load file
fileName = '191209_FTD_slice1_GroupB_1_cSpikes_L0.mat';
load(fileName);

%% get downsampled spike mat
    %get sampling frequency
    if ~exist('fs')
        fs=25000;
    else
    end

    %create spike matrix variable
    try
        spikeMatrix = full(cSpikes);
    catch
        spikeMatrix = full(mSpikes);
    end

%remove spikes from ref
spikeMatrix(:,find(channels==15))=zeros(size(spikeMatrix(:,find(channels==15))));

%downsample the spike matrix to sum spikes into 1 s time bins

    %this loop is a correction for if the length of the recording in s
    %is not a whole number; it then cuts off the additional fraction of a
    %second above the previous whole number
    if  floor(length(spikeMatrix)/fs)~=length(spikeMatrix)/fs; 
        %calculate number of samples to subtract to make 
        %length of rec in s a whole number
        n2del = fs*(length(spikeMatrix)/fs - floor(length(spikeMatrix)/fs));
        spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del-1),:);
    else
    end

recordDuration = round(length(spikeMatrix)); %in samples
downFactor = fs; %this should lead to there being 1 sample per second (1 Hz)
downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration/downFactor); 
new_fs = fs/downFactor;

%% find network spikes
%{

Aim here is to sum spikes across all channels into one vector across time
of the sum of spikes at each time point across all channels within 1s time
bins;

then take the maximum in this vector - this is the widest network spike;
can set a minimum of 3 channels to call it network spikes; if less than 3
say no network spike.

Can edit the downsampling above to adjust the bin sizes. 1 s bins considers
spikes occuring within 1 s of each other to be synchronous. 

One issue is that making these calculations with a sliding time window may
be better as spikes could occur closely in time but be in separate bins;
this is less of an issue the larger the time bins (1 s is very large in
the context of spikes). So could improve this by downsampling to 1000 Hz
then summing 50 ms time bins using sliding windows, then finding the max
no. of spikes in any given 50 ms across the array.

However this initial 1 s method works well for following up the raster
plots which use 1 s time bins

%}

%sum all electrodes across time points:
netw_spikes_vec = sum(downSpikeMatrix,2); 

down_sampled_time_of_largest_netw_spike = find(netw_spikes_vec == max(netw_spikes_vec));
%take the first time point if there are multiple time points with the
%maximum num spikes:
down_sampled_time_of_largest_netw_spike = down_sampled_time_of_largest_netw_spike(1); 

%calculate time of this network spike in recording to know TIME to plot
time_to_plot = down_sampled_time_of_largest_netw_spike * downFactor;
timeRange = time_to_plot - fs  :   time_to_plot; %to plot 1 s
%also calculate which channels participate in order to know ELECs to plot
channels_to_plot = find(downSpikeMatrix(down_sampled_time_of_largest_netw_spike,:) == 1);

%get recording file without spikes
rec_fileName = strcat(fileName(1:strfind(fileName,'Spikes')-3),'.mat');
%{
above line works by finding the position of the S in Spikes, then remove
from -3 characters including this S this which will be the m or c depending on spike
detection method and also removes the preceding underscore (e.g. removing _cS in _cSpikes).
Then at the end i concatenate '.mat')
%}

%% now plot the filtered data

load(rec_fileName);
if      strcmp(fileName(strfind(fileName,'Spikes')-1),'c')
        method = 'cwt';
        L = fileName(strfind(fileName,'Spikes')+8:end-4);
        multiplier = 500; % not needed
elseif  strcmp(fileName(strfind(fileName,'Spikes')-1),'m')
        method = 'Manuel';
        multiplier = fileName(strfind(fileName,'Spikes')+7:end-4);
        L = 500; %not needed
else
        disp('error!!!!!!!')
end


[spikeTrain, finalData, threshold] = detectSpikes(dat(timeRange,electrodes_to_plot), method, multiplier,L);

%now plot the 'finalData'...