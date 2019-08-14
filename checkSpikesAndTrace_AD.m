files = dir('*mSpikes_5.mat')

for file=1:length(files)
%% plot spike counts
%load('MPT070119_4B_DIV28_cSpikes_L0.mat')
%load('MPT070119_4B_DIV28.mat','channels')

load(files(file).name,'mSpikes','channels');
spikeMatrix=full(mSpikes);
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

sumsp=sum(adj_spikeMatrix);
count_mat=zeros(8);
count_mat([1 8 57 64])=nan;
count_mat([2:7])=sumsp(1:6);
count_mat([9:56])=sumsp(7:54);
count_mat([58:63])=sumsp(55:60);

figure
heatmap(count_mat)
title(strcat(files(file).name,{' mSpikes-5SD'}))
clear spikeMatrix
clear adj_spikeMatrix

end

%% plot grid traces to confirm channel order is correct
load('MPT070119_4B_DIV28.mat','dat')
                %filter     % currently running out of memory    
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat)); 

tracepltOrder=filteredMatrix(:,pltOrder);
option='order'; %set option to order so it plots in order of input data
%if data is not sorted into ID order, use option 'id' and input channels
%variable. Or, input data in plot order, then set option to 'id' and 
%set channels variable input as channel IDs in plot order
fs=25000;
downFactor = (fs/1000) %down factor means new n samples
%is old n samples / down factor 
% fs/1000 gives 1khz or 1 sample per ms

%option='id';
%chans=channels(pltOrder);

figure
gridTrace_AD(tracepltOrder, downFactor,option,channels,fs);


