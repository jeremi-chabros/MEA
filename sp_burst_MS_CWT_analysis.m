clear all
%if got spikes already go to spikes folder:
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\best_cSpikes_0'

%files = dir('tc*Spikes*');  % for hpc
files = dir('MPT*Spikes*');  % for ctx
    files = files(~contains({files.name}, 'adjM'));%remove unwanted files
    files = files(~contains({files.name}, 'Filt'));%remove unwanted files

sampling_fr = 25000;

progressbar
%% loop through 
for file = 1:length(files)
    %% load and get basic firing
    tic
        data=load(files(file).name,'*Spikes','channels'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
        try
            cSpikes=data.cSpikes;
            disp(strcat(files(file).name,'...cSpikes loaded'))
        catch
            mSpikes=data.mSpikes;
            disp(strcat(files(file).name,'...mSpikes loaded'))
        end
        
        try 
    channels=data.channels;
        catch
            if strcmp(files(file).name(1:2),'MP')
    channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
    channels=channels.channels;
    disp('channels loaded from dat file')
            elseif strcmp(files(file).name(end-5),'L')
    channels=load(strcat(files(file).name(1:end-15),'.mat'),'channels'); %some spike mats dont have channels var saved
    channels=channels.channels;
    disp(strcat('chans from: ',files(file).name(1:end-15),'.mat'))
            else
    channels=load(strcat(files(file).name(1:end-14),'.mat'),'channels'); %some spike mats dont have channels var saved
    channels=channels.channels;
    disp(strcat('chans from: ',files(file).name(1:end-14),'.mat'))
            end
        end
        clear data
    try
        spikeMat=full(cSpikes);
        disp('spikes based on cSpikes')
    catch
        spikeMat=full(mSpikes);
        disp('spikes based on mSpikes')
    end

spikeCounts=sum(spikeMat);
ActiveSpikeCounts=spikeCounts(find(spikeCounts>9));  %spikes of only active channels ('active'= >9)

output(file).rec = files(file).name(1:end-4); %file name
output(file).spikes = sum(spikeMat);        %cell containing n. spikes for each channel

FiringRates=spikeCounts/(length(spikeMat)/sampling_fr); %firing rate in seconds
Active_FRs=FiringRates(find(spikeCounts>9)); %firing rates in Hz of active channels only

if length(ActiveSpikeCounts)<4 
output(file).mean_FR = 0; %if there are fewer than 4 active channels, set to 0 and exclude culture
output(file).median_FR = 0;
else
output(file).mean_FR = mean(Active_FRs);
output(file).median_FR = median(Active_FRs);
end
output(file).N_active_Es = length(ActiveSpikeCounts);

%get rid of NaNs where there are no spikes; change to 0
if isnan(output(file).mean_FR);
    output(file).mean_FR=0;
else
end
if isnan(output(file).median_FR);
    output(file).median_FR=0;
else
end
disp('spike stats done')

%% get burst info
samplingRate=sampling_fr;
method ='Bakkum';
%note, Set N = 30 (min number of bursts)
%ensure bursts are excluded if fewer than 3 channels (see inside burst Detect
%function)
%to change min channels change line 207 of burstDetect.m
%to change N (min n spikes) see line 170 of burstDetect.m

[burstMatrix, burstTimes, burstChannels] = burstDetect(spikeMat, method, samplingRate);
nBursts=length(burstMatrix);
%trainCombine=sum(spikeMatrix, 2);
if length(burstMatrix)>0
    for Bst=1:length(burstMatrix)
    sp_in_bst(Bst)=sum(sum(burstMatrix{Bst,1}));
    end
    sp_in_bst=sum(sp_in_bst);
else
    disp('no bursts detected')
    sp_in_bst=0;
end

%CV of IBI
    %CV = st dev / mean
    %get IBIs
    end_times = burstTimes(1:end-1,2); %-1 because no IBI after end of last burst
    sta_times = burstTimes(2:end,1); %start from burst start time 2
    IBIs      = sta_times -end_times;
    % calculate CV of IBI and non need to convert from samples to seconds
    % (as relative measure it would be the same)
output(file).CVIBI = round((std(IBIs)/mean(IBIs)),3); %3 decimal places
    % clear unneeded vars
    clear end_times sta_times IBIs
output(file).BurstRate=60*(nBursts/(length(spikeMat(:,1))/samplingRate));
output(file).frac_in_burst = sp_in_bst/sum(sum(spikeMat));


%need to go intro burst detect and edit as it is not deleting the bursts
%with <5 channels from burstChannels and burstTimes hence they are longer
%need this for easier plotting of burst
disp('burst stats done')

clear Spikes cSpikes mSpikes
clear spikeCounts
clear ActiveSpikeCounts
clear FiringRates
clear Active_FRs
clear  active_chanIndex active_spike_mat 
clear channels   

disp('file done')
toc
progressbar(file/length(files));
end

        %% save 
        method = 'cwtANDms';
        threshold = '_params_L0AND-5SD';
        fileName = strcat(method, threshold, '_ctx_stats', '.mat'); 
        save(fileName, 'output');