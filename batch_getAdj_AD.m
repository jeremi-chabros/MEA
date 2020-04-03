% batch get adjM
clear all
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

files = dir('*FTD*Spikes*.mat*');  % where your .mat files are
files = files(~contains({files.name}, 'adjM'));
sampling_fr = 25000;
progressbar
for file = 1:length(files)
    if ~~exist(strcat(files(file).name(1:end-4), '_adjM', '.mat'))
    fprintf(strcat('\n',files(file).name(1:end-4),' \n already done','\n'))
    else
    fprintf(strcat('\n',files(file).name(1:end-4),' \n calculating adj. mat. ...','\n'))
    tic;
    % load 
        data=load(files(file).name,'*Spikes','channels'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
        
        try
    cSpikes=data.cSpikes;
        catch
    mSpikes=data.mSpikes;
        end
        
        try 
    channels=data.channels;
        catch
    channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
    channels=channels.channels;
    disp('channels loaded from dat file')
        end
        
        clear data
        try
    spikeMat=full(cSpikes);
        catch
    spikeMat=full(mSpikes);
        end

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
end