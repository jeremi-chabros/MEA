% batch get adjM

cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\best_cSpikes_0'

files = dir('*Spikes*.mat');  % where your .mat files are

sampling_fr = 25000;
progressbar
for file = 1:length(files)
    
    tic;
    % load 
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

%get adjM
method = 'tileCoef';
adjM = getAdjM(spikeMat, method, 0,0.05); %0 means no downsampling; 0.05 is sync window in s

        % save 
        fileName = strcat(files(file).name(1:end-4), '_adjM', '.mat'); 
        % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
        save(fileName, 'adjM','channels');
        toc
        progressbar(file/length(files));
        
end