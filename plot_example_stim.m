%{
example plot of a stimulation 
%}

fileName = '{190830_slice1stim5.mat}';
load(fileName);
yGap = 100; % vertical gap bewteen traces 
%electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
%    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
electrodesToPlot = [find(channels==72),find(channels==73),find(channels==74),...
    find(channels==75),find(channels==76),find(channels==77)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
%timeRange = 1: fs * 0.01;
%timeRange = 1:length(dat);
timeRange = 1298750:1299250;

figure 
for electrode = 1:length(electrodesToPlot)
            try 
            plot(filteredMatrix(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1),...
                'Color',[0,0,0])
            hold on 
            catch
            plot(filteredData(timeRange, electrodesToPlot(electrode)) - yGap * (electrode -1))
            hold on   
            end
end

aesthetics 
removeAxis 
sb=scalebar;
sb.Position=[1,min(filteredMatrix(timeRange))-(100*length(electrodesToPlot))];
sb_hoz = [int2str(sb.XLen/fs*1000),' ms']; %in ms
sb_ver = [int2str(sb.YLen),' ?V']; %in uV
sb.hTextX_Pos= [-100,-100]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-100,-100];
if strfind(fileName,'_') %remove underscores for title
    fileName(strfind(fileName,'_'))='';
else
end
title({[int2str(round(length(timeRange)/fs*1000)),'{ s of recording from }',fileName],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});