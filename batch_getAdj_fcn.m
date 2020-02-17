function batch_getAdj_fcn(method,files,sync_win,num_samples,ds_rate)
% batch get adjM
% clear all
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

%files = dir('*FTD*Spikes*.mat*');  % where your .mat files are
%files = files(~contains({files.name}, 'adjM'));
sampling_fr = 25000;
progressbar
for file = 1:length(files)
    if ~~exist(strcat(files(file).name(1:end-4), '_adjM_',num2str(sync_win/ds_rate), '.mat'))
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
    %remove spikes from ref
    spikeMat(:,find(channels==15))=zeros(size(spikeMat(:,find(channels==15))));

    % get adjM
    % method = 'tileCoef';
    if ~exist('sync_win')
        disp('need to determine sync_win!!'); %i.e. 50 ms synchroncity window is the default
    else
    end
    adjM = getAdjM(spikeMat, method, num_samples,sync_win); %0 means no downsampling; 0.05 is sync window in s
        % save 
        fileName = strcat(files(file).name(1:end-4), '_adjM_',num2str(sync_win/ds_rate), '.mat'); 
        % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
        sync_win_adjusted = sync_win/ds_rate;
        save(fileName, 'adjM','channels','sync_win_adjusted');
        toc
        clear cSpikes mSpikes channels
        progressbar(file/length(files));
    end
end
end
