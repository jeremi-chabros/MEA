function grd_elec_organoid_fcn(filename,channelIDs)

%{
This script will ground a single channel from a single file
specify these below
will work for a raw, filtered or spike data file (any spike detection
method)
%}
%clear all
%filename='SMPT190923_2B_DIV21_cSpikes_L0.1254.mat';
%channelID=46; %channel that you want to ground (use MCS coordinate)
load(filename);

if ~~exist('dat')
    
    try
        dat(:,find(channels==channelID))=dat(:,find(channels==15));%assumes 15 is the MCS ID of the reference electrode
        save(filename, 'ADCz','dat','channels','fs','header','uV','-v7.3');
    catch
        filteredMatrix(:,find(channels==channelID))=filteredMatrix(:,find(channels==15));%assumes 15 is the MCS ID of the reference electrode
        save(filename, 'channels','fs','filteredMatrix','-v7.3');
    end
    
else
    
    try
        Spikes=full(cSpikes);
    catch
        Spikes=full(mSpikes);
    end
    
    %remove spikes from ref
    Spikes(:,find(channels==15)) = zeros(length(Spikes),1);
    %set channels to = ref (i.e. 0 spikes)
    for j = 1:length(channelIDs)
        Spikes(:,find(channels==channelIDs(j)))=Spikes(:,find(channels==15));
    end
    
    %save the spike file
    if ~~exist('cSpikes')
        cSpikes=sparse(Spikes);
        save(filename, 'channels','cSpikes','thresholds','-v7.3');
        clear cSpikes
    elseif ~~exist('mSpikes')
        mSpikes=sparse(Spikes);
        save(filename, 'channels','mSpikes','thresholds','-v7.3');
        clear mSpikes
    else
        disp('error identifying spike method')
    end
    
end

end


