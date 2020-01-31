% Adjust the paramters for the Bakkum 2014 Algorithm (and the minChannel) 
% and see what happens 

% also try out different spikes; mSpikes, tSpikes, pSpikes 

% 1209 6A DIV 22

%load('/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/feature_extraction/matlab/data/goodSpikes/KO_12_09_17-6A_DIV22_info.mat')
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\cSpikes.2507'
files=dir('*.mat');
filename = files(86).name;
load(filename)
spikeMatrix = cSpikes;

%% set some paramters 
samplingRate = 25000; 



%% Look at number of detected bursts vs numChannel 

numChannel = 1:10; 
N = 10; 

numBurst = zeros(size(numChannel));
for nC = numChannel
    burstMatrix = bakkumBurstDetect(spikeMatrix, samplingRate, N, nC);
    numBurst(nC) = length(burstMatrix);
end 

figure
plot(numChannel, numBurst) 
xlabel('Minimum number of channels') 
ylabel('Number of bursts detected')
title([filename(1:9),'-',filename(11:12),'-',filename(14:18),' ','(minimum number of spikes - N - fixed to 10)'])
lineThickness(3) 
aesthetics

%% Look at number of detected bursts vs numSpike 

numSpike = 10:30; 
numBurst = zeros(size(numSpike));
nC = 1;
for nS = 1:length(numSpike)
    burstMatrix = bakkumBurstDetect(spikeMatrix, samplingRate, numSpike(nS), nC);
    numBurst(nS) = length(burstMatrix);
end 

figure
plot(numSpike, numBurst) 
xlabel('Minimum number of spikes') 
ylabel('Number of bursts detected')
title([filename(1:9),'-',filename(11:12),'-',filename(14:18),' ','(minimum number of channels - N - fixed to 1)'])
lineThickness(3) 
aesthetics

%% Look at both 
figure
%numChannel = 5:10; 
%numSpike = 10:30; 
numChannel = [3 4 5 6 7]; 
numSpike = [20 30 40 50 60];
numBurst = zeros(size(numSpike));

%run for one recording (could add to run for multiple then average then
%plot)
for nC = numChannel
    for nS = 1:length(numSpike)
        burstMatrix = bakkumBurstDetect(spikeMat, samplingRate, numSpike(nS), nC);
        numBurst(nS) = length(burstMatrix);
    end 
    plot(numSpike, 60*(numBurst/(length(spikeMat(:,1))/samplingRate)));
    xlabel('Minimum number of spikes') 
    ylabel('Burst Rate (/min)')
    lineThickness(1.5) 
    aesthetics
    hold on 
end 

lg = legend('3', '4', '5', '6', '7'); 
title(lg,'Minimum of channels')
lg.FontSize = 14;
legend boxoff
xticks(numSpike)

%% fract in bursts
figure
for nC = numChannel
    for nS = 1:length(numSpike)
        burstMatrix = bakkumBurstDetect(spikeMat, samplingRate, numSpike(nS), nC);
        numBurst(nS) = length(burstMatrix);
            for Bst=1:length(burstMatrix)
                sp_in_bst(Bst)=sum(sum(burstMatrix{Bst,1}));
            end
            sp_in_bst=sum(sp_in_bst);
            frac_in_burst(nS) = sp_in_bst/sum(sum(spikeMat));
    end 
    plot(numSpike, 100*frac_in_burst);
    xlabel('Minimum number of spikes') 
    ylabel('% Spikes within bursts')
    lineThickness(1.5) 
    aesthetics
    hold on 
end 

lg = legend('3', '4', '5', '6', '7'); 
title(lg,'Minimum of channels')
lg.FontSize = 14;
legend boxoff
xticks(numSpike)
