
%% whole raster
if ~exist('fs')
    fs=25000;
else
end

try
    spikeMatrix = full(cSpikes);
catch
    spikeMatrix = full(mSpikes);
end
%remove spikes from ref
spikeMatrix(:,find(channels==15))=zeros(size(spikeMatrix(:,find(channels==15))));
%sum(spikeMatrix(:,find(channels==15))); %check no spikes
%correction if length rec in s is not a whole number 
if ~isinteger(length(spikeMatrix)/fs) 
    %calculate number of samples to subtract to make 
    %length of rec in s a whole number
    n2del = fs*(length(spikeMatrix)/fs - round(length(spikeMatrix)/fs));
    spikeMatrix=spikeMatrix(1:length(spikeMatrix)-n2del,:);
else
end

recordDuration = length(spikeMatrix); %in samples
downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs); 

figure
imagesc(downSpikeMatrix')
%newmat=downSpikeMatrix(121:180,1:60);imagesc(newmat')%used for plotting
%only last part of spike matrix
aesthetics 
ylabel('Electrode') 
xlabel('Time (s)')
cb = colorbar;
% ylabel(cb, 'Spike count')
ylabel(cb, 'Spike Frequency (Hz)') 
cb.TickDirection = 'out';
% cb.Ticks = 0:5; % for slice 5 specifically
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)
caxis([0 5]);

yticks([1, 10:10:60])

%% log 
figure
imagesc(log(downSpikeMatrix'))

aesthetics 
ylabel('Electrode') 
xlabel('Time (s)')
cb = colorbar;
% ylabel(cb, 'Spike count')
ylabel(cb, 'Log Spike Frequency (Hz)') 
cb.TickDirection = 'out';
% cb.Ticks = 0:5; % for slice 5 specifically
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)


yticks([1, 10:10:60])
