function batchGetSpike_function(datadir,files,method,multiplier,L) 
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

    %files = dir('190830_slice1stim5.mat');  % where your .mat files are 
    %files = files(~contains({files.name}, 'Spikes'));
    %files = files(~contains({files.name}, 'Filt'));%remove unwanted files
    %disp(strcat('number of files to do = ',num2str(length(files)))
    % variable name containing MEA voltage recording should be either of 
    % these two:
    voltageVar = 'electrodeMatrix'; 
    voltageVar2 = 'dat';
    % assume it takes the form numSamp x numChannels
    % samplingRate = 25000; 
    %progressbar
    for file = 1:length(files)
        try 
            data = load(files(file).name, voltageVar); 
            data = data.(voltageVar);
            channels = load(files(file).name, 'channels'); 
            channels = channels.('channels');
        catch 
            data = load(files(file).name, voltageVar2);
            data = data.(voltageVar2);
            channels = load(files(file).name, 'channels'); 
            channels = channels.('channels');
            fprintf('Data loaded successfully \n') 
        end 
        % data = data.(voltageVar); % since matlab load struct 
        % data = electrodeMatrix
        % detect spikes
        % tic; 
        % mSpikes = sparse(getSpikeMatrix(data, 'Manuel', 5));
        % tSpikes = sparse(getSpikeMatrixAlex(data, 'Tim', 14)); %AD edited to use my getSpikeMatrix.m edited script; threshold changed from 8 to 12
        % pSpikes = sparse(getSpikeMatrix(data, 'Prez', 4));
         %L=-0.1254; %changed L to 0 to confirm it changes spike rates - it did indeed increase them as it should
         % L = log(false detection / false omission) / 36.7368
         % e.g. if you want low sensitivity (cost of false detection 1000x
         % greater), do log(1000)/36.7368 to derive L 
         
         %% get spikes and save
         if strcmp(method,'cwt')
         cSpikes = sparse(getSpikeMatrixAlex(data, 'cwt', multiplier, L)); %AD added, need to adjust loss parameter
         %toc         
         fileName = strcat(files(file).name(1:end-4), '_cSpikes_L',num2str(L), '.mat'); 
         % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
         save(fileName, 'cSpikes','channels');
         progressbar(file/length(files),[],[]);
         
         elseif strcmp(method,'Manuel')
         mSpikes = sparse(getSpikeMatrixAlex(data, 'Manuel', multiplier, L)); %AD added, need to adjust loss parameter
         %toc 
         fileName = strcat(files(file).name(1:end-4), '_mSpikes_',num2str(multiplier), '.mat'); 
         % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
         save(fileName, 'mSpikes','channels');
         progressbar(file/length(files),[],[]);
         
         elseif strcmp(method,'abs')
         aSpikes = sparse(getSpikeMatrixAlex(data, 'abs', multiplier, L));
         fileName = strcat(files(file).name(1:end-4), '_aSpikes_',num2str(multiplier), '.mat');
         save(fileName, 'aSpikes','channels');
         progressbar(file/length(files),[],[]);
         else
             disp('method inputted incorrectly')
         end 
    %D:\MECP2_2019_AD\Scripts_and_Output\S2.1.SpikeMatrix_Scripts
    %this is where current script is located and called functions are
end 