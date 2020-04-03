function batchGetSpike_function(datadir,files)
% read all .mat files in directory, look for spikeMatrix or dat
% then extract spikes, save as SPARSE MATRIX


% note that this overwrites rather than append to the info file

% I initially wanted to do this together with batch processing,
% but since I am still experimenting on the results that are obtained
% from the various spike detection algorithms, I want to save them as
% sparse spike file first, so that in the future I can play around with
% them without having to load the raw data again (unless I want to try
% out another spike detection parameter / method)



%% some parameters
%progressbar
for file = 1:length(files)
    if ~~exist(strcat(files(file).name(1:end-4),'_TTX','.mat'))
            data = load(files(file).name, 'dat');
            data = data.('dat');
            channels = load(files(file).name, 'channels');
            channels = channels.('channels');
            fprintf('Data loaded successfully \n')       
% get thresholds       
            for elec = 1:length(data,2)
            thresholds( = getSpikeMatrixTTX(data, 'Manuel', multiplier, L); %AD added, need to adjust loss parameter
            mSpikes = sparse(mSpikes);
            %toc
            fileName = strcat(files(file).name(1:end-4), '_mSpikes_',num2str(multiplier), '.mat');
            % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
            save(fileName, 'mSpikes','channels');
            progressbar(file/length(files),[],[]);
            
        %D:\MECP2_2019_AD\Scripts_and_Output\S2.1.SpikeMatrix_Scripts
        %this is where current script is located and called functions are
    else
        disp('no TTX file found')
    end
end