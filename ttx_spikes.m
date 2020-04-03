% get TTX spikes using absolute thresholds from corresponding electrodes

clear all; close all

ttxfiles = dir('*FTD*TTX.mat');                   % where your .mat files are
ttxfiles = ttxfiles(~contains({ttxfiles.name}, 'Spikes'));
ttxfiles = ttxfiles(~contains({ttxfiles.name}, 'ttx4m'));
ttxfiles = ttxfiles(~contains({ttxfiles.name}, '191210'));

refPeriod_ms = 1.5;
method = 'abs';
L = 0;
original_m = 3.5; %multiplier used for mSpikes before TTX

progressbar('files','elecs')

for ttxfile = 1:length(ttxfiles)
    %strcat(ttxfiles(ttxfile).name(1:end-4), '_cSpikes_L',num2str(L), '.mat')
    ttxfileName = ttxfiles(ttxfile).name;
    
    disp('getting thresholds from pre-ttx file')
    try
        spikefile = strcat(ttxfileName(1:end-8),'_mSpikes_',num2str(original_m),'.mat');
        load(spikefile,'thresholds');
    catch
        
        %manual correction for second ttx file
        if ~~contains(ttxfileName,'TTX2')
            spikefile = strcat(ttxfileName(1:end-10),'1_mSpikes_',num2str(original_m),'.mat');
            load(spikefile,'thresholds');
        else
            spikefile = strcat(ttxfileName(1:end-9),'1_mSpikes_',num2str(original_m),'.mat');
            load(spikefile,'thresholds');
        end
    end
    
    % loop to cancel process if file already done 
    %     if ~~exist(strcat(ttxfiles(ttxfile).name(1:end-4), '_cSpikes_L',num2str(L), '.mat'))
    if ~isempty(dir(strcat(ttxfiles(ttxfile).name(1:end-4), '_aSpikes_-*', '.mat')));
        disp('file done')
        
    else
        
        disp('loading ttx data')
        load(ttxfileName);
        
        pre_ttx_thresholds = thresholds;
        clear thresholds;
        
        disp('detecting spikes using absolute threshold')
        
            spikeMatrix = zeros(size(dat));
            finalData = zeros(size(dat)); %without initialising, you get memory error
        for elec = 1:length(pre_ttx_thresholds)
            
            multiplier = pre_ttx_thresholds(elec);
            
            % prevent spikes being detected in ref or grounded elecs
            max_spike_amplitude = -7;
            pre_ttx_thresholds(find(pre_ttx_thresholds > max_spike_amplitude)) = -100;
            
            [spikeMatrix(:, elec), finalData(:, elec), thresholds(elec)] = detectSpikes(dat(:, elec), method, multiplier,L,refPeriod_ms);
            
            %remove spikes from grounded electrodes
            spikeMatrix(:, find(pre_ttx_thresholds > max_spike_amplitude)) = 0;
            spikeMatrix(:, find(channels == 15)) = 0;
            progressbar(ttxfile/length(ttxfiles),elec/length(pre_ttx_thresholds))
            
        end
        
        % get spikes and save
        disp('saving ttx spike file')
        if strcmp(method,'cwt')
            cSpikes = sparse(spikeMatrix);
            %toc
            fileName = strcat(ttxfiles(ttxfile).name(1:end-4), '_cSpikes_L',num2str(L), '.mat');
            % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
            save(fileName, 'cSpikes','channels','thresholds');
            
        elseif strcmp(method,'Manuel')
            mSpikes = sparse(spikeMatrix);
            %toc
            fileName = strcat(ttxfiles(ttxfile).name(1:end-4), '_mSpikes_',num2str(multiplier), '.mat');
            % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
            save(fileName, 'mSpikes','channels','thresholds');
            
        elseif strcmp(method,'abs')
            aSpikes = sparse(spikeMatrix);
            %toc
            fileName = strcat(ttxfiles(ttxfile).name(1:end-4), '_aSpikes_based_on_mSpikes_',num2str(original_m),'_mthresh_',num2str(mean(pre_ttx_thresholds)), '.mat');
            save(fileName, 'aSpikes','channels','thresholds');
        else
            disp('method inputted incorrectly')
        end
        
    end
    
    progressbar(ttxfile/length(ttxfiles),0)
    
end