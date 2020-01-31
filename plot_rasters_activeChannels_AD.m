%plot multiple raster plots

%manually load desired spike mat

%channels to plot
spikeCounts=sum(full(mSpikes));
active_chanIndex=find(spikeCounts>0);
figure
for plt=1:length(active_chanIndex)
spikeTrain = spikeMat(:,active_chanIndex(plt));
option='line';
timewindow=(60000:60000+25*500)';
if sum(spikeTrain(timewindow))>5
        subplot(length(active_chanIndex),1,plt)
    singleRastPlot(spikeTrain(timewindow), option); %note there will be error if no/few spikes!
%title(strcat('channel ID:_',int2str(active_chanIndex(plt)))) 
else
    disp('no spikes')
end
end

%currently  a bug with one channel

%% one minute (compare to coterril et al)

%choose channels to plot and update activechanindex

active_chanIndex=active_chanIndex([1:3,5:8]);

figure
for plt=1:length(active_chanIndex)
spikeTrain = spikeMat(:,active_chanIndex(plt));
option='line';
%timewindow=(60000:60000+25000*60)';
timewindow=(1:length(spikeTrain))'; %entire recording
if sum(spikeTrain(timewindow))>0 %all spike counts (increase to only include those with n spikes)
        subplot(length(active_chanIndex),1,plt)
    singleRastPlot(spikeTrain(timewindow), option); %note there will be error if no/few spikes!
%title(strcat('channel ID:_',int2str(active_chanIndex(plt)))) 
else
    disp('no spikes')
end
end

sb=scalebar
sb.XLen=25000*5; %multiply by desired s
sb.XUnit='samples (5s at 25kHz)';
%sb.XUnit='5 s';
sb.YLen=0; 
%sb.YUnit = 
sb.Position=[0,-0.5];
    sb.hTextX_Pos=[-30,-100]
    sb.hTextY_Pos=[0,-1000];
   