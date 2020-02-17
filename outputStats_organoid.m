% check spikes for given set of recs ;
% input = .mat data files

%cd to data files folder
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.1.MPT190515_mat_files'

%get spike matrices (if already done comment out)
%batchGetSpike
clear all; close all
%if got spikes already go to spikes folder:
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

files = dir('*FTD*mSpikes_3.mat');  % where your .mat files are
files = files(~contains({files.name}, 'ttx'));
files = files(~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'adjM'));
files = files(~contains({files.name}, '191210'));
files = files(~contains({files.name}, 'stim'));
files = files(~contains({files.name}, 'cleaned'));

sampling_fr = 25000;

spikes_only = 0;
burst_stats = 0;
cor_ctrl = 0; % 1 means correlate randomly shuffled trains; be sure to change the sync window for the ctrl
fc_figures = 0;%change to 1 to add plots
g_figures = 0;%graph theory
g_measures = 1;

progressbar
%% loop through
for file = 1:length(files)
    %% load and get basic firing
    tic
    data=load(files(file).name,'*Spikes','channels','thresholds'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
    try
        mSpikes=data.mSpikes;
        disp('loaded mSpikes')
    catch
        cSpikes=data.cSpikes;
        disp('loaded cSpikes')
    end
    
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
    %remove ref channel spikes:
    spikeCounts(find(channels == 15)) = 0;        %cell containing n. spikes for each channel
    active_chanIndex=find(spikeCounts>=10);
    ActiveSpikeCounts=spikeCounts(active_chanIndex);  %spikes of only active channels ('active'= >9)
    
    output(file).rec = files(file).name(1:end-4); %file name
    
    %find group info
    if ~isempty(strfind(files(file).name(1:end-4),'Grp'))
        output(file).grp = output(file).rec(strfind(files(file).name(1:end-4),'Grp')+3);
    elseif ~isempty(strfind(files(file).name(1:end-4),'Group'))
        output(file).grp = output(file).rec(strfind(files(file).name(1:end-4),'Group')+5);
    else
        output(file).grp = [];
    end
    
    if   isempty(strfind(output(file).rec,'ttx')) & isempty(strfind(output(file).rec,'TTX'))
        
        output(file).ttx = 0; %1 means ttx %0 means no ttx
        
    elseif   ~isempty(strfind(output(file).rec,'ttx')) | ~isempty(strfind(output(file).rec,'TTX'))
        
        output(file).ttx = 1;
        
    else
        disp('error finding whether file involved ttx')
    end
    
    
    output(file).spikes = spikeCounts;        %cell containing n. spikes for each channel
    
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
        if burst_stats == 1;
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
        end
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
            %clear spikeMat
            clear cSpikes mSpikes
            clear burstMatrix burstChannels burstTimes Bst nBursts sp_in_bst
            
            if cor_ctrl == 1
                
                
                %                 bin_size = 0.01;
                %                 length(spikeMat(:,1))/25;
                for k = 1:100 %for running multiple randomisation iterations
                    disp(strcat(num2str(101-k),'_iterations_remaining'))
                    tic
                    for elec = 1:length(spikeMat(1,:))
                        %{
 approach:
                take each 1 ms window and randomly shuffle them
                get a vector of indices for each 1 ms period
                randomly shuffle this
                then resort the spike train using this index
                        %}
                        
                        
                        syn_win = 0.175;
                        down_factor = 25;
                        fs = 25000;
                        fr = sum(spikeMat(:,elec))/(length(spikeMat(:,elec))/fs); %firing rate %divide by seconds
                        tSim = length(spikeMat(:,elec))/fs; %duration of simulation in s
                        nTrials = 1;%number of trials
                        dt = (1/fs)*down_factor; %time bins; set to 1 ms currently
                        [rand_spike_vec, tVec] = poissonSpikeGen(fr, tSim, nTrials,dt);
                        rand_spikeMat(elec,:) = abs(rand_spike_vec);
                    end
                    rand_spikeMat = rand_spikeMat';
                    
                    %plot to compare spike mat spike counts and rand spike
                    %counts over x iterations
                    a(k,:) = sum(spikeMat) - sum(rand_spikeMat);
                    
                    adjM2 = getAdjM(rand_spikeMat, method, 0,0.007); %0 means no downsampling; 0.05 is sync window in s
                    adjM2 = adjM2(active_chanIndex,active_chanIndex);
                    adjM2 = adjM2 - eye(size(adjM2));
                    adjM2(find(isnan(adjM2)))=0;
                    
                    randSTTC(k) = round(sum(sum(adjM2))/(length(adjM2)*(length(adjM2)-1)),3);
                    
                    fileNameSpikes = strcat(files(file).name(1:end-4), '_CTRL_spikeMat.mat');
                    fileNameadjM = strcat(files(file).name(1:end-4), '_CTRL_adjM_0.175.mat');
                    save(fileNameSpikes, 'rand_spikeMat','-v7.3');
                    save(fileNameadjM, 'adjM2','-v7.3');
                    
                    clear adjM2 rand_spikeMat
                    toc
                    
                end %for running multiple randomisation iterations
                
                output(file).STTC_RCtrl= mean(randSTTC);
                output(file).nrmlsd_cnnctvty = output(file).meanSTTC/output(file).STTC_RCtrl;
                %figure; plot(mean(a))
                
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
            end
            disp('fc stats done')
            %% Graph metrics
            if g_measures == 1
                disp('beginning graph measures')
                % Basic: calculate over three absolute thresholds
                % 0.3 0.5 0.7
                count1 = 1; %to track threshold iterations
                for cutoff = [0.3 0.5 0.7]
                    
                    threshold = cutoff;
                    
                    % % choose either whole mat or only active channels
                    %   edges=adjM1;
                    edges=active_adjM;
                    
                    edges(edges < threshold) = 0.0001;
                    
                    buEdges=edges;
                    buEdges(find(buEdges==0.0001))=0;
                    buEdges(find(buEdges>0.0001))=1;
                    %remove nodes with no edges in binary network
                    active_nodes_vec=find(sum(buEdges>0));
                    buEdges=buEdges(active_nodes_vec,active_nodes_vec);
                    ActiveChanLabels=channels(active_chanIndex);
                    NodeLabels=ActiveChanLabels(active_nodes_vec);
                    
                    if length(buEdges)==0 %if there are no nodes with edges above threshold
                        netw_density(count1) = 0;%else it would be NaN with below formula
                    else
                        netw_density(count1) = sum(sum(buEdges))/(length(buEdges)*(length(buEdges)-1));
                    end
                    meanDegree(count1)       = mean(sum(buEdges));
                    netw_size (count1)       = length(buEdges); %num. nodes
                    count1 = count1 + 1;
                end
                
                output(file).netw_density   =   mean(netw_density);
                output(file).meanDegree     =   mean(meanDegree);
                output(file).netw_size      =   mean(netw_size);
                
                disp('simple graph measures done, starting complex measures')
                % Complex: calculate over three proportion thresholds
                % 40, 60 and 80 %
                count2 = 1;
                for cutoff = [40 60 80]
                    
                    threshold = prctile(adjM(:), cutoff);
                    
                    % %                     choose either whole mat or only active channels
                    %                     edges=adjM1;
                    edges=active_adjM;
                    
                    edges(edges < threshold) = 0.0001;
                    
                    buEdges=edges;
                    buEdges(find(buEdges==0.0001))=0;
                    buEdges(find(buEdges>0.0001))=1;
                    %remove nodes with no edges in binary network
                    active_nodes_vec=find(sum(buEdges>0));
                    buEdges=buEdges(active_nodes_vec,active_nodes_vec);
                    ActiveChanLabels=channels(active_chanIndex);
                    NodeLabels=ActiveChanLabels(active_nodes_vec);
                    
                    
                    
                    if round(max(sum(buEdges))/5)*5>0 % check for enough nodes
                        % clustering
                        C = clustering_coef_bu(buEdges);
                        C = mean(C);
                        % path length
                        D  = distance_bin(buEdges);
                        PL = charpath(D,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively
                        % small world
                        
                        % rich club and betweenness centrality
                        k               =   max(sum(buEdges));
                        [Rc,Nk,Ek]      =   rich_club_bu(buEdges,k);
                        RC              =   max(Rc);
                        maxKrand        =   min(find(Rc==RC));
                        
                        RC_nodes_vec    =   find(sum(buEdges) == maxKrand);
                        
                        BC_vec          =   betweenness_bin(buEdges);
                        BC              =   mean(BC_vec(RC_nodes_vec));
                        BC_norm         =   BC / ((length(BC_vec)-1) * (length(BC_vec)-2));
                        
                        % normalisation
                        disp('computing normalisations')
                        for r_i = 1:100
                            
                            [R,eff]     = randmio_und(buEdges, 100); %preserves density and distribution
                            
                            CrVEC(:,r_i)= clustering_coef_bu(R);
                            Cr(r_i)     = mean(CrVEC(:,r_i));
                            
                            Dr          = distance_bin(R);
                            PLr(r_i)    = charpath(Dr,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively
                            
                            [Rcr,Nk,Ek] = rich_club_bu(R,maxKrand);
                            RCr(r_i)    = max(Rcr);
                            
                            RCr_nodes_vec  =   find(sum(R) == maxKrand);
                            
                            BCr_vec          =   betweenness_bin(R);
                            BCr(r_i)         =   mean(BCr_vec(RCr_nodes_vec));
                            BCr_norm(r_i)    =   BCr(r_i)/((length(BCr_vec)-1) * (length(BCr_vec)-2));
                            clear BCr_vec
                            
                        end
                        disp('normalisation computed')
                        all_CC(count2)     =   C/mean(Cr);
                        all_PL(count2)     =   PL/mean(PLr);
                        all_SW(count2)     =   all_CC(count2)/all_PL(count2);
                        all_RC(count2)     =   RC/mean(RCr);
                        all_BC(count2)     =   BC/mean(BCr);
                        all_BC_n(count2)   =   BC_norm/mean(BCr_norm);
                        
                    else
                        
                        all_CC(count2)     =    0;
                        all_PL(count2)     =    length(active_chanIndex)-1;
                        all_SW(count2)     =    0;
                        all_RC(count2)     =    0;
                        all_BC(count2)     =    0;
                        all_BC_n(count2)   =    0;
                        
                    end
                    
                    count2 = count2 + 1;
                    disp('a prctle threshold iteration completed')
                    clear CrVEC 
                    
                end
                
                output(file).CC     =   mean(all_CC);
                output(file).PL     =   mean(all_PL);
                output(file).SW     =   mean(all_SW);
                output(file).RC     =   mean(all_RC);
                output(file).BC     =   mean(all_BC);
                output(file).BC_n   =   mean(all_BC_n);
                output(file).CCmax  =   all_CC(3);
                output(file).PLmax  =   all_PL(3);
                output(file).SWmax  =   all_SW(3);
                output(file).RCmax  =   all_RC(3);
                output(file).BCmax  =   all_BC(3);
                output(file).BC_nmax=   all_BC_n(3);
                
                disp('complex graph measures complete')
                
            end
            
        end %end of loop if no active channels
    end
    disp('file done')
    toc
    progressbar(file/length(files));
    
    clearvars -except files file spikes_only burst_stats cor_ctrl ...
        fc_figures g_figures g_measures sampling_fr output
    
end


%% save
disp('saving...')
method = 'manuel';
threshold = '3';
fileName = strcat(method,'_' ,threshold, '_organoid_graph_stats', '.mat');
save(fileName, 'output');
xldata = struct2table(output);
xldata(:,'spikes') = [];
writetable(xldata, strcat(fileName(1:end-4),'.xlsx'))
