% load spike matrix and compute adjacency using STTC and save for multiple
% recordings

%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\best_cSpikes_0'
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\spikes_temp'
clear all
datadir = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
cd(datadir);
files = dir('*MPT*DIV28*Spikes*het*');
%files = files(~contains({files.name}, 'het'));%remove unwanted files
files = files(~contains({files.name}, 'adjM'));%remove unwanted files


progressbar
for file=1:length(files)
    tic
    if~exist(strcat(files(file).name(1:end-4), '_adjM_DA60', '.mat'));
    load(files(file).name,'mSpikes');
    fprintf('\n \n loaded spikes \n \n \n')
        load(strcat(files(file).name(1:end-8),'.mat'),'channels');
            fprintf('\n \n loaded coords \n \n \n')
    spikeMat = full(mSpikes);
    
    %get active spike positions
    %spikeCounts=sum(spikeMat);
    %active_chanIndex=find(spikeCounts>9);
    
    %get adjM
method = 'tileCoef';
fprintf('\n \n calculating adjM \n \n \n')
%adjM = getAdjM(spikeMat(:,active_chanIndex), method, 0,0.05); %0 means no downsampling; 0.05 is sync window in s
adjM = getAdjM(spikeMat, method, 0,0.5); 
fprintf('\n \n done adjM; now saving... \n \n \n')

%save
fileName = strcat(files(file).name(1:end-4), '_adjM_DA60', '.mat'); 
        % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
        save(fileName, 'adjM','channels');
    toc
        progressbar(file/length(files));
    else
        disp('done already')
    end
end