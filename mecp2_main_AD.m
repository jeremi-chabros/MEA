%% The Mecp2 Project 

% this is the main script for the mecp2 project and aims to be a brief
% walkthrough of the usual steps taken for MEA data anlaysis. 
% Note that batch analysis is not included here yet. 
%Currently this is focused on generating some visualisations and getting an overview of the activity level of MEAs.

% a lot of this is taken directly from the script I wrote for the organoid
% project 

% Author: Tim Sit 
% Last update: 20180627

%% Add dependencies to path

% TODO: rebrand those code into some sort of MEA processing code
%addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/analysis_functions_ts/'))
% addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/organoid_data_analysis/'))

% scale bar 
%addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/chenxinfeng4-scalebar-4ca920b/'))
% heatmaps 
%addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/heatMap/'))
% human colours 
%addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/XKCD_RGB/'))
% cwt spike detection 
%addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/continuous_wavlet_transform/'))

% orgaoind project 
% addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project')); 

% figure2eps
%addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figure2epsV1-3'));

%% Convert Data from .raw to .mat 

% I assume this is done for now 

%% Perform spike detection (whole MEA implementation) 

%load a .mat data file first

method = 'Manuel';
% method = 'cwt';
multiplier = 7.5;
L = 0; 
timeRange = 1: fs * 720;
%timeRange = 900:3000; %in samples (for burst plot)
% timeRange = 110 * fs: 185 * fs -1;


%cd 'D:\MECP2_2019_AD\Scripts_and_Output\S2.1.SpikeMatrix_Scripts'
%[spikeMatrix, filteredMatrix] = getSpikeMatrixAlex(dat(timeRange, :), method, multiplier, L);

%load spike matrix
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mSpikes7.5'
%select matrix or...
load('MPT070119_1A_DIV14_mSpikes_7.5.mat','mSpikes');
spikeMatrix = full(mSpikes);


%load filtered signal matrix
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\filt_mats'
%select filtered signal; should contain fs and channels
%load('MPT070119_1A_DIV14_Filtd.mat',VAR)
% if no channels, get this from mat file load('MPT070119_1A_DIV14.mat','channels')

%% Filtered traces / grid traces (whole MEA implementation)

downFactor = 4; % down sample factor for making grid trace
cd 'D:\MECP2_2019_AD\Scripts_and_Output'
%filteredMatrix=filteredMatrix(timeRange,:);
gridTrace_AD(filteredMatrix(:,:), downFactor,'id',channels,fs) %need to add folders to path containing aesthetics function etc.


%% use below for plotting burst across network; find burst in burstmatrix
%get filteredMat (load from filt mats folder or just get)
%load dat file first
                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat)); 
%set time range to burst range and plot several plots to find a good one
%adjust length of plots to ~2 ms
%e.g. for mpt190403-6c-div14 
%timeRange = 900:3000; %in samples (for burst plot)
%set downfactor to 1 or 2
burst_to_plot =119; %look in burst channels and find one with most channels
start_samps=burstTimes(burst_to_plot,1);
%end_samps=burstTimes(burst_to_plot,2);
end_samps=start_samps+25*10 %plot time from start, set multiplier to desired ms 
timeRange =start_samps:end_samps;
downFactor = 4; % down sample factor for making grid trace
figure
gridTrace_AD(filteredMatrix(start_samps:end_samps,:), downFactor,'id',channels,fs) %need to add folders to path containing aesthetics function etc.
bcs=burstChannels{119,1};
bcsIDs=channels(bcs)';

%below plots multiple plots for a burst to find best spikes
%for plot = 1:length(timeRange)
%    figure
%start=round((plot-1)*(length(timeRange))+1);
%endp=round((plot-1)*(length(timeRange))+1+(length(timeRange)));
%gridTrace_AD(filteredMatrix(start:endp,:), downFactor,'id',channels,fs) %need to add folders to path containing aesthetics function etc.
%gridTrace_AD(filteredData, downFactor,'id',channels,fs)
%end

%% MEA Raster Plot 

figure
recordDuration = length(spikeMatrix) / fs;
downSpikeMatrix = downSampleSum(spikeMatrix, recordDuration * 1/5); 

h = imagesc(downSpikeMatrix' ./5);  %note if one channel has way more spikes then this will dominate scale bar
%change scalebar to log spikes
%h = imagesc(log(downSpikeMatrix') ./5); 
%logFR_matrix = log(downSpikeMatrix);

aesthetics 
ylabel('Electrode') 
xlabel('Time (s)')
cb = colorbar;
% ylabel(cb, 'Spike count')
ylabel(cb, 'Spike Frequency (Hz)') 
cb.TickDirection = 'out';
% cb.Ticks = 0:5; % for slice 5 specifically
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)


set(h, 'AlphaData', ~isnan(downSpikeMatrix')) % for NaN values

timeBins = 5; % 5 second separation between marks
% timePoints = 1:timeBins:floor(length(spikeMatrix) / fs); 
timePoints = 0:20:floor(length(spikeMatrix) / fs); 
yticks([1, 10:10:60])
xticks(timePoints); 
xticklabels(string(timePoints * 5));
% xticklabels(string(timePoints -1 ));

yLength = 800; 
xLength = yLength * 2; 
set(gcf, 'Position', [100 100 xLength yLength])


%% Multiple single trace plots 

yGap = 100; % vertical gap bewteen traces 
electrodesToPlot = [3, 16, 53]; % list of electrodes to plot
timeRange = 1: fs * 0.01;

figure 

for electrode = 1:length(electrodesToPlot)
            try 
            plot(filteredMatrix(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1))
            hold on 
            catch
            plot(filteredData(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1))
            hold on   
            end
end

aesthetics 
removeAxis 

%% raw data Multiple single trace plots with channels and time window and file in title - AD

yGap = 100; % vertical gap bewteen traces 
electrodesToPlot = [15, 4]; % list of electrodes to plot
time_s = 0.1 %time length to plot is s
timeRange = 1: fs * time_s; %last number here = time wish to plot in seconds
%timeRange = 1:length(dat); %last number here = time wish to plot in seconds

figure 

    data=dat;

for electrode = 1:length(electrodesToPlot)
    %x_time_ms=[1:length(timeRange)]/(fs/1000);
    x_time_s=[1:length(timeRange)]/(fs);
    y= data(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1);
            plot(x_time_s,y,'Color','b')
            hold on 

end

filename=files(2).name;

aesthetics
if length(timeRange)<25000
    plotTime=int2str(length(timeRange)/25);
    tUnit='{ms }';
else
    plotTime=int2str(length(timeRange)/25000);
    tUnit='{s }';
end

electrodeIDs = channels(electrodesToPlot)';  %may need to load channels variable
%load(strcat(files(30).name(1:end-4)),'channels'); %change file accordingly

title([plotTime,tUnit,'{of raw signal from }',filename(1:9),'-',filename(11:12),'-',filename(14:18),...
    newline,'{Electrodes IDs (top to bottom): }',int2str(electrodeIDs)]);
%use curly brackets for a string (necessary to impose finger space)
%add scalebar?? (cant change y lim due to imposed y gap)
aesthetics 
removeAxis 
       sb = scalebar;
       %sb.YLen = round(round(max(max(electrodeMatrix))-min(min(electrodeMatrix))),-1); %round to nearest 10
       %sb.XLen = length(electrodeMatrix)/fs; 
       sb.YUnit = '\muV';
       %sb.XUnit = 'ms';
       sb.XUnit = 's';
       sb.Position = [-time_s/10, min(y)-40];
       %may need to adjust accoridng to legnth of plot window
       %sb.hTextX_Pos = [length(trace)/10,-sb.YLen/5];
       %sb.hTextY_Pos = [-length(trace)/6,sb.YLen/10];
       
       
%% add a point to indicate the spike time
stimes = find(s(timeRange,4)==1); %search for nearest spike within plot window
% add one, may need to do manually
x_secs = stimes/fs;
y_spikes = ones(size(x))*min(y)-yGap/8;
plot(x_secs,y_spikes,'^','LineWidth',5,'MarkerFaceColor','r','MarkerEdgeColor','r');


 
%%  Multiple single trace plots filtered data

yGap = 100; % vertical gap bewteen traces 
electrodesToPlot = [15, 49, 4]; % list of electrodes to plot
timeRange = 1: fs * 10; %last number here = time wish to plot in seconds

figure 

if ~exist('filteredMatrix')
    filteredMatrix=filteredData;
else
end

for electrode = 1:length(electrodesToPlot)
            plot(filteredMatrix(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1))
            hold on 

end

filename=files(2).name;

aesthetics
if length(timeRange)<25000
    plotTime=int2str(length(timeRange)/25);
    tUnit='{ms }';
else
    plotTime=int2str(length(timeRange)/25000);
    tUnit='{s }';
end

electrodeIDs = channels(electrodesToPlot)';  %may need to load channels variable
load(strcat(files(30).name(1:end-4),'_FiltD'),'channels'); %change file accordingly

title([plotTime,tUnit,'{of filtered signal from }',filename(1:9),'-',filename(11:12),'-',filename(14:18),...
    newline,'{Electrodes IDs (top to bottom): }',int2str(electrodeIDs)]);
%use curly brackets for a string (necessary to impose finger space)
%add scalebar?? (cant change y lim due to imposed y gap)
aesthetics 
removeAxis 

%% MEA Spike Sum Heat Map 
    %note tim's heat map is in wrong positions - needs adjustment to put
    %channels in correct place

%get spike mat
%spike_dir = 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mSpikes7.5'
%files=dir('*mSpikes*');
%spike_dir = 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\cSpikes.2507'
spike_dir = 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\cSpikes.0';
cd(spike_dir);
files=dir('*cSpikes*');

for file=1:length(files)
filename=files(file).name;
load(filename);
%spikeMatrix=full(mSpikes);
spikeMatrix=full(cSpikes);

    %load channel variable
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mats';
load(strcat(filename(1:18),'.mat'),'channels')

%adjust spike matrix from being in channel order (see channels file '1 2 3 4 etc')
%into  order of channel IDs counting along rows (by column) i.e. '21 31 41
%51 61 71 81 12 22 32 etc.'

%fast way:
pltOrder=[find(channels==21),find(channels==31),find(channels==41),... %count across columns (subplot index plots across columns, e.g. sublot (4,4,2) the plots in column 2 not row 2
        find(channels==51),find(channels==61),find(channels==71),find(channels==12),...
        find(channels==22),find(channels==32),find(channels==42),find(channels==52),...
        find(channels==62),find(channels==72),find(channels==82),find(channels==13),...
        find(channels==23),find(channels==33),find(channels==43),find(channels==53),...
        find(channels==63),find(channels==73),find(channels==83),find(channels==14),...
        find(channels==24),find(channels==34),find(channels==44),find(channels==54),...
        find(channels==64),find(channels==74),find(channels==84),find(channels==15),...
        find(channels==25),find(channels==35),find(channels==45),find(channels==55),...
        find(channels==65),find(channels==75),find(channels==85),find(channels==16),...
        find(channels==26),find(channels==36),find(channels==46),find(channels==56),...
        find(channels==66),find(channels==76),find(channels==86),find(channels==17),...
        find(channels==27),find(channels==37),find(channels==47),find(channels==57),...
        find(channels==67),find(channels==77),find(channels==87),find(channels==28),...
        find(channels==38),find(channels==48),find(channels==58),find(channels==68),...
        find(channels==78)];
adj_spikeMatrix=spikeMatrix(:,pltOrder);

%slow way takes twice as long
%output_channel_order=   [21,31,41,51,61,71,...
%                        12,22,32,42,52,62,72,82,...
%                        13,23,33,43,53,63,73,83,...
%                        14,24,34,44,54,64,74,84,...
%                        15,25,35,45,55,65,75,85,...
%                        16,26,36,46,56,66,76,86,...
%                       17,27,37,47,57,67,77,87,...
%                        28,38,48,58,68,78];

%for iteration=1:length(spikeMatrix(1,:));
%    currentChannelID=channels(iteration);
%    output_matrix_position=find(currentChannelID==output_channel_order);
%    adj_spikeMatrix(:,output_matrix_position)=spikeMatrix(:,iteration); 
%    %disp({[int2str(length(spikeMatrix(1,:))-iteration),' electrodes remaining']})
%end

figure
makeHeatMap(adj_spikeMatrix, 'rate') %choose 'rate' or 'count'
set(gcf, 'Position', [100, 100, 800, 800 * 1])
title([filename(1:9),'-',filename(11:12),'-',filename(14:18)]);
%need to set scalebar limit to 5Hz - make universal
%caxis([0,5])%comment out when checking cultures; turn on for thesis figure

clear 'adj_spikeMatrix';
cd(spike_dir);
clear 'channels';
disp({[int2str(length(files)-file),' files remaining']})

end