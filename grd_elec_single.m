%{
This script will ground a single channel from a single file
specify these below
will work for a raw, filtered or spike data file (any spike detection
method)
%}
clear all
filename='SMPT190923_2B_DIV21_cSpikes_L0.1254.mat';
channelID=46; %channel that you want to ground (use MCS coordinate)
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
        cSpikes=full(cSpikes);
        cSpikes(:,find(channels==channelID))=cSpikes(:,find(channels==15));%assumes 15 is the MCS ID of the reference electrode
        cSpikes=sparse(cSpikes);
        save(filename, 'channels','cSpikes','-v7.3');
    catch
        mSpikes=full(mSpikes);
        mSpikes(:,find(channels==channelID))=mSpikes(:,find(channels==15));%assumes 15 is the MCS ID of the reference electrode
        mSpikes=sparse(mSpikes);
        save(filename, 'channels','mSpikes','-v7.3');
    end
    
end


