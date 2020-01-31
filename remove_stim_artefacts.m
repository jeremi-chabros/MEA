% remove stim artefacts
clear all
%% identify stim artefacts and count them
fileName = '191209_FTD_slice9_GroupB_3_stim.mat';
load(fileName);
electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
%electrodesToPlot = [find(channels==12),find(channels==13),find(channels==14),...
%    find(channels==16),find(channels==17)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
%timeRange = 1: fs * 0.01;
timeRange = 1:length(dat);
%timeRange = 75084:75584;
yGap = 200; % vertical gap bewteen traces 

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
sb_ver = [int2str(sb.YLen),' uV']; %in uV
sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-1000,-1000];
    fileName1=fileName;
if strfind(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))='';
    fileName1=strcat('{',fileName1,'}');
    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
else
    title({[int2str(round(length(timeRange)/fs*1000)),'{ ms of recording from }',fileName],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
end


%% get the stim artefact times from channel with the number of stim artefacts identified from above

load('191209_FTD_slice9_GroupB_3_stim_mSpikes_6.mat')
sum(full(mSpikes));
electrode_to_plot=42; %not MCS ID but which column the desired E is in
multiplier=6;
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'Manuel', multiplier,0);
sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=40;
%n_spikes_to_plot=length(sp_times);
for i=1:n_spikes_to_plot
sp_peak_time=find(finalData(sp_times(i):sp_times(i)+25)==min(finalData(sp_times(i):sp_times(i)+25)));%note changed min to max as with filtered data spike is +ve
%plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25));
plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25),...
'Color',[0.5 0.5 0.5],'LineWidth',3); %all grey
hold on
ave_trace(:,i)=finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25);
end
ave_trace=(sum(ave_trace'))/(length(ave_trace(1,:)));
%plot(ave_trace,'LineWidth',8,'Color','w');
%plot(ave_trace,'LineWidth',3,'Color','r');
plot(ave_trace,'LineWidth',3,'Color','k'); %black line instead
hold off
aesthetics
box off
%change axes to voltage normalised
%change x axis to time rather than samples
xticks(linspace(0,50,3));
xticklabels(linspace(0,50,3)/25);
xlabel('time (ms)');
ylim([-50 50])
yticks(linspace(-80,80,(160/20)+1));
ylabel('filtered signal (\muV)');
set(gca,'fontsize',16)
threshvec=ones(length(ave_trace))*threshold;
hold on
plot(1:length(ave_trace),threshvec,'LineStyle','--','Color','b','LineWidth',3)
hold off
%calculate the average and overlay that onto the figure as a red line
%surrounded by white
clear ave_trace
%sp_times
stim_times_s=[round([sp_times/25000],0)]'; %stim times in s
%% remove stim artefacts from v trace
clearvars -except fileName sp_times
load(fileName);
% remove +/- 100 ms either side of spike times
ms_to_remove = 10;
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
    
%% plot to show stims removed
clear all
fileName = '191209_FTD_slice9_GroupB_3_stim_cleaned.mat';

load(fileName);
electrodesToPlot = [find(channels==82),find(channels==83),find(channels==84),...
    find(channels==85),find(channels==86),find(channels==87)]; % list of electrodes to plot
%electrodesToPlot = [find(channels==12),find(channels==13),find(channels==14),...
%    find(channels==16),find(channels==17)]; % list of electrodes to plot

                fs=25000;
                lowpass = 600; 
                highpass = 8000; 
                wn = [lowpass highpass] / (fs / 2); 
                filterOrder = 3; %changed from 3 to 5 by AD; and back
                [b, a] = butter(filterOrder, wn); 
                filteredMatrix = filtfilt(b, a, double(dat));
                
%timeRange = 1: fs * 0.01;
timeRange = 1:length(dat);
%timeRange = 75084:75584;
yGap = 200; % vertical gap bewteen traces 

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
sb_ver = [int2str(sb.YLen),' uV']; %in uV
sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
sb.hTextY_Pos= [-1000,-1000];
    fileName1=fileName;
if strfind(fileName,'_') %remove underscores for title
    fileName1(strfind(fileName1,'_'))='';
    fileName1=strcat('{',fileName1,'}');
    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName1],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
else
    title({[int2str(round(length(timeRange)/fs)),'{ s of recording from }',fileName],...
    ['{ (scalebars: horizontal }',sb_hoz,'{ vertical }',sb_ver,')']});
end
