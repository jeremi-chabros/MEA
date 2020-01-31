% remove stim artefacts

%{ 

Last edit: Alex Dunn 17-12-2019

Department of Phys., Dev. and Neuroscience,
University of Cambridge

This script will remove artefacts from MEA recordings on an MEA2100 MCS
system and save the output with the original file name plus 'cleaned'.
This will also save images showing the removal of artefacts - please check
these to confirm it has worked.

This assumes mSpikes_6 files have already been created (i.e. spike
detection using the 'Manuel' method and a multiplier of 6. File names need
to be saved as the original file names with the voltage trace plus
'_mSpikes_6.mat'

Note this script assumes that the manipulation of voltage lasts no more
than around 3 ms per injection of current (or voltage deflection)

ensure that you set out at the top which files you want to clean

%}

%% determine stim files; 
clear all;close all; clc
files = dir('*FTD*Group*stim*.mat*');             
files = files(~contains({files.name}, 'Spikes'));
files = files(~contains({files.name}, 'cleaned'));

progressbar('files','sections')

for file = 1:length(files)  
    % identify file to do
    fileName = files(file).name;

if ~exist(strcat(fileName(1:end-4),'_cleaned.mat'))
    %% identify stim artefacts and count them
    load(fileName);
    electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
    %electrodesToPlot = [find(channels==12),find(channels==13),find(channels==14),...
    %find(channels==16),find(channels==17)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
    %timeRange = 1: fs * 0.01;
    timeRange = 1:length(dat);
    %timeRange = 74959:75709;
    yGap = 200; % vertical gap between traces 

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
sb.Position=[1,min(filteredMatrix(timeRange))-(yGap*length(electrodesToPlot))];
%sb_hoz = [int2str(sb.XLen/fs*1000),' ms']; %in ms
sb_hoz = [int2str(sb.XLen/fs),' s']; %in s
sb_ver = [int2str(sb.YLen),' uV']; %in uV
sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-1000,-1000];
    fileName1=fileName;
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName1=strcat('{',fileName1(1:end-4),'}');
%    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName1],...
%    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
else
%    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName],...
%    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
end

    %make title font smaller
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.7;

    %save as PNG
    fprintf(strcat('\n','\n',files(file).name(1:end-4),' saving fig...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_whole_rec.png'));    
    close all; clear cSpikes mSpikes 

progressbar([],1/6)
    
    %% plot a 30 ms close up of the stim artefacts 
    load(fileName);
    electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
    %electrodesToPlot = [find(channels==12),find(channels==13),find(channels==14),...
    %find(channels==16),find(channels==17)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
    % find first spike 
    load(strcat(fileName(1:end-4),'_mSpikes_6.mat'))
    sumspikes = sum(full(mSpikes));
    modespikes = mode(sumspikes);
    chans_with_mode = find(sumspikes(electrodesToPlot)==modespikes); %looks which electrode to plot out of the ones plotted above only
    electrode_to_plot = electrodesToPlot(chans_with_mode(1));
    sp_times_plot = find(mSpikes(:,electrode_to_plot)==1);
    
    if isempty(sp_times_plot)
        fprintf(strcat(' \n \n \n \n \n \n \n \n \n \n','\t \t <strong> warning: </strong> \n  \n \t cannot find stim artefacts in \n','\t',...
            fileName,'\n \t need to try template method \n \t or manually remove spikes based on timing... \n \n \n \n \n \n \n \n \n \n'))
    else
        timeRange_samples = 750/2; %750/2 will mean plotting 15 ms before and after if its 25 samples per ms
        
        %added correction so that timerange cannot be negative i.e. if
        %first spike is very early in recording
        if  sp_times_plot(1) > timeRange_samples
            timeRange = sp_times_plot(1) - timeRange_samples : sp_times_plot(1) + timeRange_samples;
        else
            sp_times_plot = sp_times_plot(find(sp_times_plot>timeRange_samples))
            timeRange = sp_times_plot(1) - timeRange_samples : sp_times_plot(1) + timeRange_samples;
        end
            
        yGap = 200; % vertical gap between traces 

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
        sb.Position=[1,min(filteredMatrix(timeRange))-(yGap*length(electrodesToPlot))];
        sb_hoz = [int2str(sb.XLen/fs*1000),' ms']; %in ms
        %sb_hoz = [int2str(sb.XLen/fs),' s']; %in s
        sb_ver = [int2str(sb.YLen),' uV']; %in uV
        sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
        sb.hTextY_Pos= [-1000,-1000];
        fileName1=fileName;
        
        if  contains(fileName,'_') %remove underscores for title
            fileName1(strfind(fileName1,'_'))=' ';
            fileName1=strcat('{',fileName1(1:end-4),'}');
            title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName1],...
            ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
            %    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
            %    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
        else
            title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName],...
            ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
            %    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
            %    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
        end

        %make title font smaller
        ax = gca;
        ax.TitleFontSizeMultiplier = 0.7;

        %save as PNG
        fprintf(strcat('\n','\n',files(file).name(1:end-4),' saving fig...', '\n','\n'))
        saveas(gcf,strcat(fileName(1:end-4),'_30msCloseUp.png'));    
        close all; clear cSpikes mSpikes 
    

    progressbar([],2/6)

%% get the stim artefact times from channel with the number of stim artefacts identified from above
    % and plot the overlay of spikes in a channel with stim-like spikes
close all;
load(strcat(fileName(1:end-4),'_mSpikes_6.mat'))
sumspikes = sum(full(mSpikes));
modespikes = mode(sumspikes);

chans_with_mode = find(sumspikes(electrodesToPlot)==modespikes); %looks which electrode to plot out of the ones plotted above only
electrode_to_plot = electrodesToPlot(chans_with_mode(1));
if length(chans_with_mode) == 0
    chans_with_mode = find(sumspikes==modespikes);               %looks at all electrodes
    electrode_to_plot = electrodesToPlot(chans_with_mode(1));
    disp(strcat({'plotted channel: '}, int2str(channels(1))))
else
    disp(strcat({'plotted channel: '}, int2str(electrodesToPlot(1))))
end
    
% plot the spikes in that electrode
multiplier=6;
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'Manuel', multiplier,0);

%%%                         %%%           
%%%                         %%%
%%%                         %%%
sp_times=find(spikeTrain==1);
%%%                         %%%
%%%                         %%%
%%%                         %%%

figure

if sum(spikeTrain)>=40
    n_spikes_to_plot=40;
else
    n_spikes_to_plot=sum(spikeTrain);
end

%n_spikes_to_plot=length(sp_times);
for i=1:n_spikes_to_plot
    %AD 17/12/19: changed to find negative peak up to 50 (rather than 25) samples (2 ms rather than 1 ms) 
    % from spike time this is because stim can last longer - this script wont work for all stim durations
    
sp_peak_time=find(finalData(sp_times(i):sp_times(i)+50)==min(finalData(sp_times(i):sp_times(i)+50)));%note changed min to max as with filtered data spike is +ve
%plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25));
plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25),...
'Color',[0.5 0.5 0.5],'LineWidth',3); %all grey
hold on
all_trace(:,i)=finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25);
end
ave_trace=(sum(all_trace'))/(length(all_trace(1,:)));
%plot(ave_trace,'LineWidth',8,'Color','w');
%plot(ave_trace,'LineWidth',3,'Color','r');
plot(ave_trace,'LineWidth',3,'Color','k'); %black line instead
hold off
aesthetics
box off
%change axes to voltage normalised
%change x axis to time rather than samples
xticks(linspace(0,50,3));
xticklabels([linspace(0,50,3)/25]-1);
xlabel('time relative to negative peak (ms)');

%calibrate y axis
nearestValue = 100; % i.e. nearest mutliple of 100
ymax = ceil(max(max(all_trace))/nearestValue)*nearestValue;
ymin = ceil(min(min(all_trace))/nearestValue)*nearestValue;
if abs(ymax) >= abs(ymin)
    yval = abs(ymax);
else
    yval = abs(ymin);
end

ylim([-yval yval])

if yval >= 300
    increment = 100;
else
    increment = 50;
end

yticks(linspace(-yval,yval,((2*yval)/increment)+1));
ylabel('filtered signal (\muV)');
set(gca,'fontsize',16)
threshvec=ones(length(ave_trace))*threshold;
hold on
plot(1:length(ave_trace),threshvec,'LineStyle','--','Color','b','LineWidth',3)
hold off
%calculate the average and overlay that onto the figure as a red line
%surrounded by white
clear ave_trace all_trace
%sp_times
stim_times_s=[round([sp_times/25000],0)]'; %stim times in s

%save the figure
    fileName1=fileName;
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName1=strcat(fileName1(1:end-4));
    title({[int2str(n_spikes_to_plot),' spikes overlaid from channelXY ',...
        int2str(channels(electrode_to_plot)),' in ',fileName1],...
    [' ']});
else
    title({[int2str(n_spikes_to_plot),' spikes overlaid from channelXY ',...
        int2str(channels(electrode_to_plot)),' in ',fileName1],...
    [' ']});
end

    %make title font smaller
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.7;

    %save as PNG
    fprintf(strcat('\n','\n',files(file).name(1:end-4),' saving fig...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot)),...
        '_stim_artefact_overlay.png'));    
    close all; clear cSpikes mSpikes 

progressbar([],3/6)


%% remove stim artefacts from v trace

load(fileName);
% remove +/- 100 ms either side of spike times
ms_to_remove = 100;
n_samples = ms_to_remove*fs/1000;  %number of samples to remove before and after stim time
for i=1:size(dat,2)
    for j=1:length(sp_times)
        try
            dat(sp_times(j)-n_samples:sp_times(j)+n_samples,i) = ... 
            dat(sp_times(j)-n_samples:sp_times(j)+n_samples,find(channels==15));
        catch
            if     sp_times(j)-n_samples < 1
                
                dat(1:sp_times(j)+n_samples,i) = ... 
                dat(1:sp_times(j)+n_samples,find(channels==15));
            
            elseif sp_times(j)+n_samples > length(dat)
                
                dat(sp_times(j)-n_samples:end,i) = ... 
                dat(sp_times(j)-n_samples:end,find(channels==15));
            
            else
                disp('error in code used to remove artefacts')
            end
        end
    end
    %disp('electrode done')
end
%save new file
disp('saving...')
                save(strcat(fileName(1:end-4),'_cleaned.mat'), 'ADCz','channels','dat','fs','uV', '-v7.3');    
              
progressbar([],4/6)

%% plot to show stims removed
close all
clear dat channels fs
cleaned_fileName = strcat(fileName(1:end-4),'_cleaned.mat');
    load(cleaned_fileName);
    electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
    %electrodesToPlot = [find(channels==12),find(channels==13),find(channels==14),...
    %find(channels==16),find(channels==17)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
    %timeRange = 1: fs * 0.01;
    timeRange = 1:length(dat);
    %timeRange = 74959:75709;
    yGap = 200; % vertical gap between traces 

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
sb.Position=[1,min(filteredMatrix(timeRange))-(yGap*length(electrodesToPlot))];
%sb_hoz = [int2str(sb.XLen/fs*1000),' ms']; %in ms
sb_hoz = [int2str(sb.XLen/fs),' s']; %in s
sb_ver = [int2str(sb.YLen),' uV']; %in uV
sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-1000,-1000];
    fileName1=cleaned_fileName;
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName1=strcat('{',fileName1(1:end-4),'}');
%    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName1],...
%    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
else
%    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName],...
%    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
end

    %make title font smaller
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.7;

    %save as PNG
    fprintf(strcat('\n','\n',files(file).name(1:end-4),' saving fig...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_whole_rec_after.png'));    
    close all; clear cSpikes mSpikes 

progressbar([],5/6)
    
    %% plot a 30 ms close up of the stim artefacts having been removed
    load(cleaned_fileName);
    electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
    %electrodesToPlot = [find(channels==12),find(channels==13),find(channels==14),...
    %find(channels==16),find(channels==17)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
    % find first spike 
    load(strcat(fileName(1:end-4),'_mSpikes_6.mat'))
    sumspikes = sum(full(mSpikes));
    modespikes = mode(sumspikes);
    chans_with_mode = find(sumspikes(electrodesToPlot)==modespikes); %looks which electrode to plot out of the ones plotted above only
    electrode_to_plot = electrodesToPlot(chans_with_mode(1));
    sp_times_plot = find(mSpikes(:,electrode_to_plot)==1);
    timeRange_samples = 750/2; %750/2 will mean plotting 15 ms before and after if its 25 samples per ms
    
        if  sp_times_plot(1) > timeRange_samples
            timeRange = sp_times_plot(1) - timeRange_samples : sp_times_plot(1) + timeRange_samples;
        else
            sp_times_plot = sp_times_plot(find(sp_times_plot>timeRange_samples))
            timeRange = sp_times_plot(1) - timeRange_samples : sp_times_plot(1) + timeRange_samples;
        end
    
    yGap = 200; % vertical gap between traces 

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
sb.Position=[1,min(filteredMatrix(timeRange))-(yGap*length(electrodesToPlot))];
sb_hoz = [int2str(sb.XLen/fs*1000),' ms']; %in ms
%sb_hoz = [int2str(sb.XLen/fs),' s']; %in s
sb_ver = [int2str(sb.YLen),' uV']; %in uV
sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-1000,-1000];
    fileName1=cleaned_fileName;
if  contains(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))=' ';
    fileName1=strcat('{',fileName1(1:end-4),'}');
    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
%    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
%    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
else
    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
%    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
%    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
end

    %make title font smaller
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.7;

    %save as PNG
    fprintf(strcat('\n','\n',files(file).name(1:end-4),' saving fig...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_30msCloseUpAfter.png'));    
    close all; clear cSpikes mSpikes     

    end
progressbar([],6/6)

else
    disp(strcat(fileName,'_already done, file skipped...'))
    progressbar(file/length(files),6/6);    
end

progressbar(file/length(files),[])

end

%{ 

Now check the figures.
If they look like cleaning was successful, run spike detection on the
cleaned files and see if spikes can be detected

%}