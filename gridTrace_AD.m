function gridTrace(electrodeMatrix, downFactor,option,channels,fs)

%{
load a .mat file created from mea batch convert using whole ooption

inputs

if option is 'id' then will plot channels in correct position on grid
if it's 'order' then will plot from 1-60, then need to check channel
variable to know identity of electrode

channels = vector of channel IDs in their order in the electrode
matrix]e.g. the first channel in the elec matrix is channel 47
this is due to the file conversion output

fs is smapling freq

%}


    % plot the spike trace of the 60 electrodes 
    % exect electrode matrix to be of dimension samples x nChannels
    % specify 8 x 8 grid 
    % y lims set to min/max of entire electrode matrix
    numRow = 8; 
    numColumn = 8; 
    trace = downsample(electrodeMatrix, downFactor);
    % 
    %for excuding channels... would need to edit below to account for
    %missing channels
    pL = 1:(size(electrodeMatrix, 2)+4);
    pL = pL(pL~=1); 
    pL = pL(pL~=8); 
    pL = pL(pL~=57);
    pL = pL(pL~=64);

    
    if strcmp(option, 'order') 
           for plotN = 1:size(electrodeMatrix, 2)
       subplot(numRow, numColumn, pL(plotN)) %for excluding channels
       %subplot(numRow, numColumn, plotN)
       plot(trace(:, plotN))
       ylim([min(min(trace)), max(max(trace))/2])
       title(num2str(plotN))
       axis tight
       aesthetics()
       removeAxis()
       hold on
           end 
           
    elseif strcmp(option, 'id') 
           for plotN = 1:size(electrodeMatrix, 2)
               currentChannelID=channels(plotN);
               channel_positions=   [0,21,31,41,51,61,71,0 ,...
                                    12,22,32,42,52,62,72,82,...
                                    13,23,33,43,53,63,73,83,...
                                    14,24,34,44,54,64,74,84,...
                                    15,25,35,45,55,65,75,85,...
                                    16,26,36,46,56,66,76,86,...
                                    17,27,37,47,57,67,77,87,...
                                     0,28,38,48,58,68,78,0 ];
               position_to_plot=find(channel_positions==currentChannelID); 
               %note this is column wise
               %[currentChannelID,position_to_plot]%use this print current
               %channel number and its grid identity
       subplot(numRow, numColumn, position_to_plot)
       plot(trace(:, plotN))
       axis tight %may need to take off
       ylim([min(min(electrodeMatrix)), max(max(electrodeMatrix))]) % turn off for burst plot
       %title(strcat(num2str(plotN),'/',num2str(currentChannelID))) %turn this on to check channels are plotted in correct place
       title(num2str(currentChannelID))
       aesthetics()
       removeAxis()
       hold on
           end 
    else
        display('need to choose plotting option')
    end
    
    %add scalebar

   %% Scalebar 
   
  
       subplot(numRow, numColumn, 64)
       ylim([min(min(trace)), max(max(trace))/2])
       %linkaxes %relies on some other functions from matlab
       axis off
       set(gca,'YTickLabel',[]);
       set(gca,'XTickLabel',[]);
       set(gca,'ytick',[]);
       set(gca,'xtick',[]);
       set(gca,'linewidth',3);
       sb = scalebar;
       sb.YLen = round(max(max(trace))/5,-1); %round to nearest 10
       sb.XLen = length(trace)/5; 
       sb.YUnit = '\muV';
       sb.XUnit = 's'; 
       sb.Position = [0, -0.5*sb.YLen];
       %may need to adjust accoridng to legnth of plot window
       sb.hTextX_Pos = [length(trace)/10,-sb.YLen/5];
       sb.hTextY_Pos = [-length(trace)/6,sb.YLen/10];
end 