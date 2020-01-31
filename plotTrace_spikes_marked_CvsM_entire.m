% plotting trace with spikes indicated

%load spike mat and/or voltage trace

%and/or get spikes and filtered matrix
electrode_to_plot=46; %not MCS ID but which column the desired E is in
%dat=dat(:,electrode_to_plot);  %second value here is the channel number
L=-0.1254;
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'cwt', 0,L);

%% concatonate pre and post stim then plot together
%add blank space in the middle to show time gap before stim
load('190705_organoidspine_slice2.mat')
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'cwt', 0,L);
pre_spikeTrain=spikeTrain(1:1500000,1);
pre_finalData=finalData(1:1500000,1);

load('190830_slice1stim5.mat')
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'cwt', 0,L);

%% whole recording
fs=25000;
down_sample_factor=25;
NewSampsN = length(dat(:,1))/down_sample_factor; %divide by 25 to get 1kHz
New_fs = fs/down_sample_factor;
dss=downSampleSum(finalData,NewSampsN);
dssSpikesC=downSampleSum(spikeTrain,NewSampsN);
figure
subplot(2,1,1)
plot(finalData)
aesthetics
box off
axis off
sb=scalebar
       sb.hTextX_Pos = [-1000000000 1];
       sb.hTextY_Pos = [-1000000000 1];
       sb.YLen = round(round(max(max(finalData))-min(min(finalData)))/5,-1); %',1' means round to nearest 10
       sb.XLen = length(finalData)/10; 
       sb.Position = [450000 30]
       %label scalebar
       title([int2str(length(finalData)/fs),' s from channel ',int2str(channels(electrode_to_plot)),...
           ' (scalebars: horizontal = ',int2str(sb.YLen),' uV; vertical = ', ...
           int2str(length(finalData)/10/fs),' s)'])
       


subplot(2,1,2)
singleRastPlot(spikeTrain, 'line')
%ylim([-50 50])

%% spikes closer whole rec
fs=25000;
down_sample_factor=25;
NewSampsN = length(dat(:,1))/down_sample_factor; %divide by 25 to get 1kHz
New_fs = fs/down_sample_factor;
dss=downSampleSum(finalData,NewSampsN);
dssSpikesC=downSampleSum(spikeTrain,NewSampsN);
figure
%subplot(2,1,1)
plot(finalData)
dave=ylim;
aesthetics
box off
axis off
sb=scalebar
       sb.hTextX_Pos = [-1000000000 1];
       sb.hTextY_Pos = [-1000000000 1];
       sb.YLen = round(round(max(max(finalData))-min(min(finalData)))/10); %',1' means round to nearest 10
       sb.XLen = round(length(finalData)/10/25000,-1)*25000; 
       sb.Position = [450000 -25]
%dave2=ylim;
%ylim(dave2)
%subplot(2,1,2)
%sp_tr_enlarge=spikeTrain*29;
%plot(sp_tr_enlarge,'.','LineWidth',5,'MarkerFaceColor','r','MarkerEdgeColor','k');
%singleRastPlot(sp_tr_enlarge, 'line')
spikePos = find(spikeTrain == 1); 
plot([spikePos spikePos], [29 39], 'k')
ylim(dave+[0 10]);

%% blown up

desired_time_to_plot_in_secs=0.02;
downsample='no'; %need to fix 'yes' version — include bit to plot ds trace and ds spiketrain!
%find part with most spikes
sp_tr_length=length(spikeTrain);

if strcmp(downsample,'no')
    n_timewindows=sp_tr_length/(desired_time_to_plot_in_secs*fs);
    samplewindow=desired_time_to_plot_in_secs*fs;
elseif strcmp(downsample,'yes')
    n_timewindows=sp_tr_length/(desired_time_to_plot_in_secs*New_fs);
    samplewindow=desired_time_to_plot_in_secs*New_fs;
else 
end



count=0;
for i=1:n_timewindows
    twindow=[1+count*samplewindow:samplewindow*i];
    count=count+1;
    sp_train=spikeTrain(twindow);
    spikes(:,i)=sum(sp_train);
end

maximum = find(spikes==max(spikes));
sample_window_highest_spikes = [1+((maximum-1)*samplewindow):samplewindow*maximum]';

%
%if strcmp(downsample,'no')
%timeWindow=[1:fs*desired_time_to_plot_in_secs]';
figure
%subplot(2,1,2)
plot(finalData(sample_window_highest_spikes))
aesthetics
axis off
sb=scalebar
       %sb.YLen = round(round(max(max(finalData))-min(min(finalData)))/10); %round to nearest 10
       sb.YLen = 5;
       sb.XLen = length(sample_window_highest_spikes)/10; 
       %sb.YUnit = '\muV';
       %sb.XUnit = 'ms';
       %sb.XUnit = 's';
       sb.Position = [250.0000  -20];
       %may need to adjust accoridng to legnth of plot window
       sb.hTextX_Pos = [-1000 -1000];
       sb.hTextY_Pos = [-1000 -1000];
%subplot(2,1,1)
%plot(spikeTrain(sample_window_highest_spikes))
%singleRastPlot(spikeTrain(sample_window_highest_spikes), 'line')
spikePos = find(spikeTrain(sample_window_highest_spikes) == 1); 
plot([spikePos,spikePos]', [20 25], 'k') %add ' to flip x axis matrix if only 2 spikes


%
%{
note: do not use downsamplesum function to downsample voltage trace
this function does not downsample, it sums time bins
if using this to downsample spike trains, you need to find
spikeTrain>1 and set it back to 1
%}

%elseif strcmp(downsample,'yes')
%timeWindow=[1:New_fs*desired_time_to_plot_in_secs]';
%figure
%subplot(2,1,1)
%plot(dss(timeWindow))
%aesthetics
%subplot(2,1,2)
%singleRastPlot(dssSpikesC(timeWindow), 'line')
%else
    disp('error')
%end
    

%% old code... 
sumC=sum(full(cSpikes));
sumM=sum(full(mSpikes));

dssSpikesC=downSampleSum(full(cSpikes),NewSampsN);
dssSpikesM=downSampleSum(full(mSpikes),NewSampsN);



subplot(4,4,1)
