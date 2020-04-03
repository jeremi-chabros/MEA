% check spikes for given set of recs ; 
% input = .mat data files

%cd to data files folder
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.1.MPT190515_mat_files'

%get spike matrices (if already done comment out)
%batchGetSpike
clear all
%if got spikes already go to spikes folder:
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\best_cSpikes_0'

files = dir('*Spikes_L-0.*.mat');  % where your .mat files are
    files = files(~contains({files.name}, 'tc'));%remove unwanted files %ctx only

sampling_fr = 25000;

fc_figures = 0;%change to 1 to add plots
g_figures = 0;%graph theory

progressbar
%% loop through 
for file = 1:length(files)
    %% load and get basic firing
    tic
        data=load(files(file).name,'*Spikes','channels'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
        cSpikes=data.cSpikes;
        try 
    channels=data.channels;
        catch
    channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
    channels=channels.channels;
    disp('channels loaded from dat file')
        end
        clear data

spikeMat=full(cSpikes);
spikeCounts=sum(spikeMat);
ActiveSpikeCounts=spikeCounts(find(spikeCounts>9));  %spikes of only active channels ('active'= >9)

output(file).rec = files(file).name(1:18); %file name
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

%% functional connectivity
method = 'tileCoef';
%get adjM - ideally already created by batch get adjM to save time when re
%analysing
active_chanIndex=find(spikeCounts>9);
    
if length(active_chanIndex)<3 %set all output to 0 etc. if no activity
                                    %else run through analyses
                                    
     output(file).meanSTTC=0;
     output(file).STTC_RCtrl=0;
     output(file).netw_density=0;
     output(file).meanDegree =0;
     output(file).CC=0;
     output(file).PL=58;
     output(file).SW=0;
     output(file).RC=0;
                                        
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
output(file).meanSTTC=sum(sum(active_adjM))/(length(active_adjM)*(length(active_adjM)-1));


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
    output(file).STTC_RCtrl=sum(sum(adjM2))/(length(adjM2)*(length(adjM2)-1));

    
    if fc_figures ==1
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

%% graph measures
    %do analysis on active channels active_adjM or adjM1 full adj mat?
    %latter
%take 95th percentile / get binary undirected adjM
 %threshold = prctile(adjM1(:), 90); % get the 95th percentile 
threshold = 0.5;

%edges=adjM1;
edges=active_adjM; %choose either whole mat or only active channels

edges(edges < threshold) = 0.0001; 
% plot the diagram of network

buEdges=edges;
buEdges(find(buEdges==0.0001))=0;
buEdges(find(buEdges>0.0001))=1;
%remove nodes with no edges in binary network
active_nodes_vec=find(sum(buEdges>0));
buEdges=buEdges(active_nodes_vec,active_nodes_vec);
ActiveChanLabels=channels(active_chanIndex);
NodeLabels=ActiveChanLabels(active_nodes_vec);


%calclate average node degree and density (density is proportion of edges
%that exist) of active channels; threshold the active adj matrix
%threshold = prctile(active_adjM(:), 90);
%threshold = 0.5;
     %active_edges=active_adjM;
%active_edges(active_edges < threshold) = 0.0001; 
%active_edges(find(buEdges==0.0001))=0;
%active_edges(find(buEdges>0.0001))=1;

%buEdges=buEdges(active_chanIndex,active_chanIndex); %need this if
%edges are whole mat not just active
if length(buEdges)==0 %if there are no nodes with edges above threshold
    output(file).netw_density=0;%else it would be NaN with below formula
else
output(file).netw_density=sum(sum(buEdges))/(length(buEdges)*(length(buEdges)-1));
end
output(file).meanDegree = mean(sum(buEdges));

%plot degree distribution / calculate node degree, cc, pl etc relative to n
%active channels, not all 60 (input matrix should be active channels only
if g_figures==1 & round(max(sum(buEdges))/5)*5>0
figure
hist(sum(buEdges));    %ensure active channels only
title('degree distribution')
xlabel('node degree')
ylabel('frequency')
xlim([0,round(max(sum(buEdges))/5)*5]) %rounded to nearest 5
aesthetics

%active weighted mat
figure
imagesc(active_adjM)
cb=colorbar
cb.Limits=[min(min(active_adjM)),1]
    xlabel('Channel number');
        ylabel('Channel number');
colormap 'cool'
MAP = colormap; %displays current colour map
MAP2=1-colormap;
colormap(MAP2)
aesthetics
box on
set(gca,'TickLength',[0 0])
%set(gca,'xtick',[])
xticks(linspace(1,length(buEdges),length(buEdges))) %adjust to n active channels
yticks(linspace(1,length(buEdges),length(buEdges))) %adjust to n active channels
%xticklabels(channels(active_chanIndex))
%yticklabels(channels(active_chanIndex))
xticklabels(NodeLabels)
yticklabels(NodeLabels)


%figure
%hist(sum(buEdges));         % or plot all channels
%title('degree distribution')
%xlabel('node degree')
%ylabel('frequency')
%xlim([0,round(max(sum(buEdges))/5)*5]) %rounded to nearest 5
%aesthetics

figure
imagesc(buEdges)
colormap([1 1 1; 0 0 0]);
    xlabel('Channel ID');
        ylabel('Channel ID');
aesthetics
box on
set(gca,'TickLength',[0 0])
%set(gca,'xtick',[])
xticks(linspace(1,length(buEdges),length(buEdges))) %adjust to n active channels
yticks(linspace(1,length(buEdges),length(buEdges))) %adjust to n active channels
xticklabels(NodeLabels)
yticklabels(NodeLabels)

else
    disp('graph figures turned off')
end

%% calculate cc
if round(max(sum(buEdges))/5)*5>0
 %generate random graph; divide cc by CCr   
C=clustering_coef_bu(buEdges);
%C=mean(C(active_chanIndex)); %if buEdges was set to include all channels
C=mean(C);
alpha=50; %of edges to rewire

%for randinho =1:100
%R = randomizer_bin_und(buEdges,alpha); % use this to preserve degree distribution
%HOWEVER, does not work with low no. nodes
    %if not preserved, change alpha
    
    %randomise preserving no. nodes and edges (i.e. density but not degree distribution)
%    Nnodes=length(buEdges);
%    Nedges=sum(sum(buEdges))/2; %note need to divide by 2 or random edge mat will have 2x as many connections
%R = makerandCIJ_und(Nnodes,Nedges);
%Cr=clustering_coef_bu(R);
%Cr_V(randinho)=mean(Cr);
%clear R
%lear Cr
%end
%Cr=mean(Cr_V);
%clear Cr_V

[R,eff]=randmio_und(buEdges, 100) %preserves density and distribution
Cr=clustering_coef_bu(R);
Cr=mean(Cr);

output(file).CC=C/Cr;

    if g_figures==1
    figure
    hist(sum(R))
    title('randomised degree distribution')
    xlabel('node degree')
    ylabel('frequency')
    xlim([0,round(max(sum(buEdges))/5)*5]) %rounded to nearest 5
    aesthetics

    figure
    imagesc(R)
    colormap([1 1 1;1 0.2 0.2])
    xlabel('Channel ID');
        ylabel('Channel ID');
    aesthetics
    box on
    set(gca,'TickLength',[0 0])
    %set(gca,'xtick',[])
    xticks(linspace(1,length(buEdges),length(buEdges))) %adjust to n active channels
    yticks(linspace(1,length(buEdges),length(buEdges))) %adjust to n active channels
    xticklabels(NodeLabels)
    yticklabels(NodeLabels)
    else 
    end

else
   output(file).CC=0;
end

%% calculate pl
if round(max(sum(buEdges))/5)*5>0
 %generate random graph; divide pl by PLr   
D = distance_bin(buEdges);
PL=charpath(D,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively

%random controls

%for randinho =1:100
%R = randomizer_bin_und(buEdges,alpha); % check node degree distribution
    %if not preserved, change alpha

    %randomise preserving no. nodes and edges (i.e. density but not degree distribution)
%    Nnodes=length(buEdges);
%    Nedges=sum(sum(buEdges))/2; %note need to divide by 2 or random edge mat will have 2x as many connections
%R = makerandCIJ_und(Nnodes,Nedges);
%
%D = distance_bin(R);
%PLr_vec(randinho)=charpath(D,0,0);
%clear R
%clear D
%end
%PLr=mean(PLr_vec);

Dr = distance_bin(R);
PLr=charpath(Dr,0,0); % 0 and 0 excludes self connections and infinities (non connected) respectively


output(file).PL=PL/PLr;
else
 output(file).PL=length(active_chanIndex)-1; %could instead set to the number of nodes -1 
end

%% small world

if round(max(sum(buEdges))/5)*5>0
output(file).SW=(output(file).CC)/(output(file).PL);
else
    output(file).SW=0;
end

if g_figures==1 & round(max(sum(buEdges))/5)*5>0
    figure 
    %plotAdj(adjM1,1:60); %plots not in topographic order...
    %adjust adjM to topographic order
    topo_order_vector=[find(channels==12),find(channels==13),find(channels==14),...
        find(channels==15),find(channels==16),find(channels==17),find(channels==21),...
        find(channels==22),find(channels==23),find(channels==24),find(channels==25),...
        find(channels==26),find(channels==27),find(channels==28),find(channels==31),...
        find(channels==32),find(channels==33),find(channels==34),find(channels==35),...
        find(channels==36),find(channels==37),find(channels==38),find(channels==41),...
        find(channels==42),find(channels==43),find(channels==44),find(channels==45),...
        find(channels==46),find(channels==47),find(channels==48),find(channels==51),...
        find(channels==52),find(channels==53),find(channels==54),find(channels==55),...
        find(channels==56),find(channels==57),find(channels==58),find(channels==61),...
        find(channels==62),find(channels==63),find(channels==64),find(channels==65),...
        find(channels==66),find(channels==67),find(channels==68),find(channels==71),...
        find(channels==72),find(channels==73),find(channels==74),find(channels==75),...
        find(channels==76),find(channels==77),find(channels==78),...
        find(channels==82),find(channels==83),find(channels==84),find(channels==85),...
        find(channels==86),find(channels==87)];
    
        %new order:
    %figure
    %imagesc(adjM1(topo_order_vector,topo_order_vector))
    %xticks(linspace(5,60,12))
    %yticks(linspace(5,60,12))
    %numerical_order=[12:17,21:28,31:38,41:48,51:58,61:68,71:78,82:87];
    %    %xticklabels(numerical_order(linspace(5,60,12)))
    %    %yticklabels(numerical_order(linspace(5,60,12)))
    %   %all labels
    %xticks(linspace(1,60,60))
    %yticks(linspace(1,60,60))
    %xticklabels(numerical_order(linspace(1,60,60)))
    %yticklabels(numerical_order(linspace(1,60,60)))
    %xlabel('Channel ID')
    
        %original order:
    %figure
    %imagesc(adjM1) %check it works...
    %ticks(linspace(5,60,12))
    %yticks(linspace(5,60,12))
    %    %xticklabels(channels(linspace(5,60,12)))
    %    %yticklabels(channels(linspace(5,60,12)))
    %    %all labels
    %xticks(linspace(1,60,60))
    %yticks(linspace(1,60,60))
    %xticklabels(channels(linspace(1,60,60)))
    %ticklabels(channels(linspace(1,60,60)))
    %label('Channel ID')

    %it works
    
    adjM1_topoOrder=adjM1(topo_order_vector,topo_order_vector);
    figure
    plotAdj(adjM1_topoOrder,1:60);
else
    disp('g figs turned off or no edges')
    
end

% plot a random graph degree distribution against a real degree dist.
%% rich club
if round(max(sum(buEdges))/5)*5>0

k=max(sum(buEdges));
[Rc,Nk,Ek] = rich_club_bu(buEdges,k);
RC=max(Rc);
maxKrand = min(find(Rc==RC));
%random controls

%for randinho =1:100
%R = randomizer_bin_und(buEdges,alpha); % check node degree distribution
    %if not preserved, change alpha
%randomise preserving no. nodes and edges (i.e. density but not degree distribution)
%    Nnodes=length(buEdges);
%    Nedges=sum(sum(buEdges))/2; %note need to divide by 2 or random edge mat will have 2x as many connections
%R = makerandCIJ_und(Nnodes,Nedges);
    

%[Rc,Nk,Ek] = rich_club_bu(R,maxKrand);
%RCr_vec(randinho)=max(Rc);
%end
%RCr=mean(RCr_vec);
[Rcr,Nk,Ek] = rich_club_bu(R,maxKrand);
RCr=max(Rcr);


output(file).RC=RC/RCr;
else 
output(file).RC=0;    
end
disp('graph stats done')



    end %end of loop if no active channels
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
        method = 'cwt';
        threshold = '0';
        fileName = strcat(method, threshold, '_stats_FINAL_ko', '.mat'); 
        save(fileName, 'output');
        