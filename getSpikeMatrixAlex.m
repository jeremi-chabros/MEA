function [spikeMatrix,finalData] = getSpikeMatrixAlex(data, method, multiplier,L)
% Loop through detectspikes for each channel 
% Assume input matrix: numSamp x numChannels
if ~exist('L')
    L = 0; 
end 
    
spikeMatrix = zeros(size(data)); 
finalData = zeros(size(data)); %without initialising, you get memory error
for j = 1:size(data, 2)
    [spikeMatrix(:, j), finalData(:, j), threshold] = detectSpikes(data(:, j), method, multiplier,L);
    % fprintf(num2str(j)) % this was for debugging, to see which electrode 
    % is making an error
end 

end 