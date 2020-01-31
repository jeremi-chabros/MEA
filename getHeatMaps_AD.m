files=dir('*SMPT*_2B*Spikes*_L0.1254*.mat*');

for file=1:length(files)
filename=files(file).name;
load(filename);
%spikeMatrix=full(mSpikes);
spikeMatrix=full(cSpikes);

    %load channel variable
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mats';
%load(strcat(filename(1:18),'.mat'),'channels')

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

option = 'logc';
figure
makeHeatMap_AD(adj_spikeMatrix, option) %choose 'rate' or 'count' or 'logc'
set(gcf, 'Position', [100, 100, 800, 800 * 1])
%title([filename(1:9),'-',filename(11:12),'-',filename(14:18)]);
%title(filename(1:end-15))
title(filename)
%need to set scalebar limit to 5Hz - make universal
if strcmp(option,'rate')
    caxis([0,5]);%comment out when checking cultures; turn on for thesis figure
elseif strcmp(option,'count')
    caxis([0, 2000]); %
elseif strcmp(option,'logc') %log10 count
    caxis([0,4]); %4 is 10 000 spikes
end
clear 'adj_spikeMatrix';
%cd(spike_dir);
clear 'channels';
disp({[int2str(length(files)-file),' files remaining']})

end