% check spikes for given set of recs ; 
% input = .mat data files

%cd to data files folder
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.1.MPT190515_mat_files'

%get spike matrices (if already done comment out)
%batchGetSpike
clear all
%if got spikes already go to spikes folder:
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

files = dir('*200114*mSpikes_5.mat');  % where your .mat files are
    %files = files(~contains({files.name}, 'tc'));%remove unwanted files %ctx only

sampling_fr = 25000;

spikes_only = 0;
fc_figures = 0;%change to 1 to add plots
g_figures = 0;%graph theory

progressbar
%% loop through 
for file = 1:length(files)
    %% load and get basic firing
    tic
        data=load(files(file).name,'*Spikes','channels'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
        mSpikes=data.mSpikes;
        try 
    channels=data.channels;
        catch
    channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
    channels=channels.channels;
    disp('channels loaded from dat file')
        end
        clear data

spikeMat=full(mSpikes);
spikeCounts=sum(spikeMat);
active_chanIndex=find(spikeCounts>=5);
ActiveSpikeCounts=spikeCounts(active_chanIndex);  %spikes of only active channels ('active'= >9)

                output(file).rec = files(file).name(1:end-4); %file name
                output(file).grp = output(file).rec(strfind(files(file).name(1:end-4),'Grp')+3);
                              
    if   isempty(strfind(output(file).rec,'ttx')) & isempty(strfind(output(file).rec,'TTX'))
    
                output(file).ttx = 0; %1 means ttx %0 means no ttx
                
    elseif   ~isempty(strfind(output(file).rec,'ttx')) | ~isempty(strfind(output(file).rec,'TTX'))
    
                output(file).ttx = 1; 
                
    else
                disp('error finding whether file involved ttx')
    end
    
    
output(file).spikes = sum(spikeMat);        %cell containing n. spikes for each channel

FiringRates=ActiveSpikeCounts/(length(spikeMat)/sampling_fr); %firing rate in seconds

%if length(ActiveSpikeCounts)<4 
%output(file).mean_FR = 0; %if there are fewer than 4 active channels, set to 0 and exclude culture
%output(file).median_FR = 0;
%else
output(file).mean_FR = round(mean(FiringRates),3);
output(file).sem_FR = round(std(FiringRates)/(sqrt(length(ActiveSpikeCounts))),3);
output(file).median_FR = round(median(FiringRates),3);
output(file).iqr_FR = round(iqr(FiringRates),3);
%end
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
if spikes_only ~= 1
    
    %turn warning off
    warning('off','MATLAB:nearlySingularMatrix');
    
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
output(file).BurstRate=round(60*(nBursts/(length(spikeMat(:,1))/samplingRate)),3);
output(file).frac_in_burst = round(sp_in_bst/sum(sum(spikeMat)),3);


%need to go intro burst detect and edit as it is not deleting the bursts
%with <5 channels from burstChannels and burstTimes hence they are longer
%need this for easier plotting of burst
disp('burst stats done')

%% functional connectivity
method = 'tileCoef';
%get adjM - ideally already created by batch get adjM to save time when re
%analysing
    
if length(active_chanIndex)<3 %set all output to 0 etc. if no activity
                                    %else run through analyses
                                    
     output(file).meanSTTC=0;
     output(file).STTC_RCtrl=0;
     %output(file).netw_density=0;
     %output(file).meanDegree =0;
     %output(file).CC=0;
     %output(file).PL=58;
     %output(file).SW=0;
     %output(file).RC=0;
                                        
     disp('no active channels')
   
    
else
        
        active_spike_mat = spikeMat(:,active_chanIndex);

        if ~exist(strcat(files(file).name(1:end-4), '_adjM', '.mat'))

    active_adjM=getAdjM(active_spike_mat, method, 0,0.05); 
    active_adjM=active_adjM-eye(size(active_adjM));
    active_adjM(find(isnan(active_adjM)))=0;

%adjM = getAdjM(spikeMat, method, 0,0.05); %0 means no downsampling; 0.05 is sync window in s
 

        else
    adjM=load(strcat(files(file).name(1:end-4), '_adjM', '.mat'));
    adjM=adjM.adjM;
    adjM1= adjM-eye(size(adjM));
    adjM1(find(isnan(adjM1)))=0;
    active_adjM=adjM1(active_chanIndex,active_chanIndex); %for extracting
    %part of full adjmat but now just compute active adjmat
    

        end

%mean correlation of active channels
output(file).meanSTTC= round(sum(sum(active_adjM))/(length(active_adjM)*(length(active_adjM)-1)),3);
output(file).semSTTC= round(std(nonzeros(triu(active_adjM)))/(sqrt(length(nonzeros(triu(active_adjM))))),3);

%% fc control - randomise FR of each active electrode and compute STTC
    %create randomised controls for each recording and save their control
    %adjMs using the batch get adjM
    
    %load spike mat and find active channels

    %save memory; clear some vars
    clear spikeMat
    clear cSpikes
    clear burstMatrix burstChannels burstTimes Bst nBursts sp_in_bst
    
    for tarly = 1:length(active_spike_mat(1,:))
               
        rand_spike_vec=active_spike_mat(:,tarly);
        rand_spike_vec=rand_spike_vec(randperm(length(rand_spike_vec)));
        rand_spike_mat(:,tarly)=rand_spike_vec;
        %check distribution
        %sts=find(rand_spike_vec==1); %spike times
        %sts2=sts(2:length(sts));       %spike times t+1
        %ISIsss=sts2-sts(1:end-1);          %spike time t+1 - t
        %ISIsss=ISIsss/25000;               %get the ISIs in seconds
        %figure
        %hist(ISIsss)
                            %distribution looks like poison process,
                            %positively skewed as mean ISI is low, with a
                            %few large ISIs
                               %if mean ISI was higher i.e. low mean Firing
                               %rate, poisson process would mean less
                               %positively skewed, closer to normal but
                               %still positvely skewed a bit
    end 
    adjM2 = getAdjM(rand_spike_mat, method, 0,0.05); %0 means no downsampling; 0.05 is sync window in s
    adjM2 = adjM2 - eye(size(adjM2));
    
    %mean STTC of randomised spike trains with same Fire rate
    output(file).STTC_RCtrl= round(sum(sum(adjM2))/(length(adjM2)*(length(adjM2)-1)),3);

    
    if fc_figures == 1
    %compare
    figure;
    imagesc(active_adjM,[min(min(active_adjM)),1]); %set colouring limits!
    xlabel('Channel ID');
        ylabel('Channel ID');
    aesthetics
    box off
    set(gca,'TickLength',[0 0])
    %set(gca,'xtick',[])
    xticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
    yticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
    xticklabels(channels(active_chanIndex))
    yticklabels(channels(active_chanIndex))
    cbFC=colorbar
    cbFC.Limits=[min(min(active_adjM)),1]


    figure;
    imagesc(adjM2,[min(min(active_adjM)),1])
    aesthetics
    box off
    set(gca,'TickLength',[0 0])
        xticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
    yticks(linspace(1,length(active_adjM),length(active_adjM))) %adjust to n active channels
    xticklabels(channels(active_chanIndex))
    yticklabels(channels(active_chanIndex))
    cbFD=colorbar
    cbFD.Limits=[min(min(active_adjM)),1]
    else 
    end

disp('fc stats done')

    end %end of loop if no active channels
end
disp('file done')
toc
progressbar(file/length(files));

clear Spikes
clear spikeCounts
clear ActiveSpikeCounts
clear FiringRates
clear Active_FRs
clear rand_spike_vec
clear rand_spike_mat
clear buEdges active_adjM active_chanIndex active_spike_mat adjM1 adjM2 
clear channels   edges    tarly 
end


        %% save 
        method = 'manuel';
        threshold = '5';
        fileName = strcat(method,'_' ,threshold, '_organoid_stats', '.mat'); 
        save(fileName, 'output');
        xldata = struct2table(output);
        xldata(:,'spikes') = [];
        writetable(xldata, strcat(fileName(1:end-4),'.xlsx'))
        