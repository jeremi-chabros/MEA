%% plot filtered trace of spikes overlaid and peaks aligned
function plotSpikes_fcn(fileName,method,parameter,refPeriod_ms,option)
% set chosen time windows, position 1 is in s, pos 2 and 3 are in ms
time_wins           = [60 100 30];
multiplier          = parameter;
L                   = parameter;

%% load data

load(fileName)
if ~~exist('aSpikes')
    spikeMatrix = full(aSpikes);
    
elseif ~~exist('cSpikes')
    spikeMatrix = full(cSpikes);
    
elseif ~~exist('mSpikes')
    spikeMatrix = full(mSpikes);
    
else
    disp('error - cant get spike mat')
end
clear aSpikes cSpikes mSpikes
disp('loaded spikes')

if strcmp(option,'spikiest')
    %find channel index (note this is the column number of the desired data,
    %not the MCS XY coordinate)
    %not MCS ID but which column the desired E is in
    electrode_to_plot   = find(sum(spikeMatrix) == max(sum(spikeMatrix)));
    %correction if it's a draw (i.e. >1 electrodes has max no. spikes)
    if  length(electrode_to_plot)>1
        electrode_to_plot = electrode_to_plot(1);
    end
    
elseif strcmp(option,'diagonal')
    electrode_to_plot = [find(channels == 22),find(channels == 33),find(channels == 44),...
        find(channels == 55),find(channels == 66),find(channels == 77)];
else
    disp('cannot identify plotting option')
end

for ei = 1:length(electrode_to_plot)
    %% load raw data / get filtered data to plot
    spike_suffix_index = strfind(fileName,'Spikes')-3; %-3 to remove _aS in aSpikes for example
    spike_suffix = fileName(spike_suffix_index+1:end-4);
    raw_dat_fileName = strcat(fileName(1:spike_suffix_index),'.mat');
    
    if          ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
            '_spike_overlay.png')) ...
            & ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
            '_',num2str(time_wins(1)),'s_trace_spikes_marked.png')) ...
            & ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
            '_',num2str(time_wins(2)),'ms_trace_spikes_marked.png')) ...
            & ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
            '_',num2str(time_wins(3)),'ms_trace_spikes_marked.png'))
        disp('all figs already done')
    else
        
        if  ~~exist(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'))
            load(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'));
            finalData =     filteredMatrix(:,electrode_to_plot(ei));
            clear filteredMatrix
            spikeTrain =    spikeMatrix(:,electrode_to_plot(ei));
        else
            load(raw_dat_fileName);
            % loop to calibrate absolute threshold for ttx files
            if contains(fileName,'TTX') | contains(fileName,'ttx')
                multiplier = thresholds(electrode_to_plot(ei));
            end
            [spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot(ei)),method, multiplier,L,refPeriod_ms);
        end
        disp('acquired filtered data and spike trains')
        %AD: adding loop as error caused if there are no spikes to plot
        if sum(spikeTrain) < 5
            disp('no spikes found')
        else
            
            if ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_spike_overlay.png'))
                
                fprintf(strcat('\n \n',fileName(1:end-4),'\n',...
                    '_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_spike_overlay.png','\n \n already done! \n \n'))
            else %create plot if not already done:
                
                %% get spike times and plot trace
                
                sp_times=find(spikeTrain==1);
                figure
                n_spikes_to_plot=50;
                %added correction if num spikes is fewer than desired number:
                if  sum(spikeTrain) < n_spikes_to_plot
                    n_spikes_to_plot = sum(spikeTrain);
                end
                
                for i=1:n_spikes_to_plot
                    sp_peak_time=find(finalData(sp_times(i):sp_times(i)+25)==min(finalData(sp_times(i):sp_times(i)+25)));%note changed min to max as with filtered data spike is +ve
                    %plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25));
                    
                    % added correction for if there is an early, such that spike time -
                    % minus 25 samples is negative, it will not plot
                    % or if any part of the trace to plot is above 50 or below -50, it
                    % will exclude this trace from the plot and the average
                    % calculation
                    if (sp_times(i))+sp_peak_time-25 < 1 | ...
                            ~isempty(find(...
                            finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25) > 50)) | ...
                            ~isempty(find(...
                            finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25) < -50))
                        % don't plot; don't add average, first line of array will be
                        % excluded by MATLAB automatically
                    else
                        plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25),...
                            'Color',[0.5 0.5 0.5],'LineWidth',3); %all grey
                        hold on
                        all_trace(:,i)=finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25);
                    end
                end
                
                ave_trace=(sum(all_trace'))/(length(all_trace(1,:)));
                %plot(ave_trace,'LineWidth',8,'Color','w');
                %plot(ave_trace,'LineWidth',3,'Color','r');
                plot(ave_trace,'LineWidth',3,'Color','k'); %black line instead
                hold off
                aesthetics
                
                %change axes to voltage normalised
                %change x axis to time rather than samples
                xticks(linspace(0,50,3));
                xticklabels([linspace(0,50,3)/25]-1);
                xlabel('time relative to negative peak (ms)');
                
                %calibrate y axis
                nearestValue = 50; % i.e. nearest mutliple of 50
                ymax = ceil(max(max(all_trace))/nearestValue)*nearestValue;
                ymin = ceil(min(min(all_trace))/nearestValue)*nearestValue;
                if  abs(ymax) >= abs(ymin)
                    yval = abs(ymax);
                else
                    yval = abs(ymin);
                end
                
                ylim([-yval yval])
                
                if      yval >= 300
                    increment = 100;
                elseif  yval >= 100
                    increment = 50;
                else
                    increment = 25;
                end
                
                yticks(linspace(-yval,yval,((2*yval)/increment)+1));
                ylabel('filtered signal (\muV)');
                set(gca,'fontsize',16)
                
                % plot the threshold if using threshold-based method
                if      strcmp(method,'Manuel') | strcmp(method,'abs')
                    threshvec = ones(length(ave_trace))*threshold;
                    hold on
                    plot(1:length(ave_trace),threshvec,'LineStyle','--','Color','b','LineWidth',3)
                    hold off
                    %calculate the average and overlay that onto the figure as a red line
                    %surrounded by white
                else
                end
                
                %% sort title and save the figure
                disp('saving plot...')
                
                fileName1=fileName;
                if  contains(fileName,'_') %remove underscores for title
                    fileName1(strfind(fileName1,'_'))=' ';
                    fileName1=strcat(fileName1(1:end-4));
                    title({[int2str(n_spikes_to_plot),' spikes overlaid from channel (XY): ',...
                        int2str(channels(electrode_to_plot(ei))),' in '],[fileName1],...
                        [' ']});
                else
                    title({[int2str(n_spikes_to_plot),' spikes overlaid from channel (XY): ',...
                        int2str(channels(electrode_to_plot(ei))),' in '],[fileName1],...
                        [' ']});
                end
                
                %make title font smaller
                ax = gca;
                ax.TitleFontSizeMultiplier = 0.7;
                
                %save as PNG
                fprintf(strcat('\n','\n',fileName(1:end-4),' saving fig...', '\n','\n'))
                saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_spike_overlay.png'));
                
            end
            
            
            
            %% spikiest 60 s
            close all
            
            time_win        = 	time_wins(1); % in s
            
            %check if already done
            if ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'s_trace_spikes_marked.png'))
                
                fprintf(strcat('\n \n',fileName(1:end-4),'\n',...
                    '_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'s_trace_spikes_marked.png \n \n already done! \n \n'))
                
            else %if not already done, carry out below...
                %if we don't have the spike train because above figure was already
                %complete, get it now...
                if  ~exist('spikeTrain')
                    spike_suffix_index = strfind(fileName,'Spikes')-3; %-3 to remove _aS in aSpikes for example
                    spike_suffix = fileName(spike_suffix_index+1:end-4);
                    raw_dat_fileName = strcat(fileName(1:spike_suffix_index),'.mat');
                    
                    if  ~~exist(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'))
                        load(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'));
                        finalData =     filteredMatrix(:,electrode_to_plot(ei));
                        clear filteredMatrix
                        spikeTrain =    spikeMatrix(:,electrode_to_plot(ei));
                        threshold = threshold;
                    else
                        load(raw_dat_fileName);
                        [spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot(ei)),method, multiplier,L,refPeriod_ms);
                    end
                    disp('acquired filtered data')
                end
                
                
                if ~exist('fs') %correction if fs isn't loaded properly (sampling freq)
                    fs = 25000;
                end
                sample_win      = time_win * (fs); % because sampling freq, fs, is in seconds
                slidage         = 10 * fs; % fs = 1 s; so this slides every 10 s
                num_windows     = (length(spikeTrain)- sample_win) / slidage;
                num_done        = 0;
                
                %progressbar
                for i = 1:num_windows %-1 so that i can start from 0
                    spike_sums(i)   = sum(spikeTrain(1+slidage*num_done:slidage*num_done+sample_win));
                    num_done        = num_done + 1;
                    %progressbar(i/num_windows)
                end
                most_spikes_index = find(spike_sums == max(spike_sums))-1; % -1 is necessary
                %because if the first window had the most spikes, we need to start from
                %one; if the 2nd window had the most spikes, we need to start from
                %1+slidage*1
                
                if  length(most_spikes_index) > 1
                    most_spikes_index = most_spikes_index(2); %if it's a draw just take second one
                else
                end
                
                %get index of samples to plot
                samples_index = 1+slidage*most_spikes_index:slidage*most_spikes_index+sample_win;
                
                %plot
                figure
                x1 = (1:length(samples_index)) / (fs); %time vec in s
                y1 = finalData(samples_index);
                plot(x1,y1,'Color','k');
                aesthetics
                box off
                %ylim([min(y1) max(y1)])
                ylim([-50 50])
                hold on
                %y2 = mean(y1)+6*std(y1);
                y2 = 50;
                sp_times = find(spikeTrain(samples_index) == 1);
                plot(x1(sp_times),y2,'v','LineWidth',0.5,'MarkerFaceColor','r','MarkerEdgeColor','r');
                
                removeAxis
                sb=scalebar;
                sb.Position=[1,-40];
                sb.XLen = 5; % in s
                sb.YLen = 5; % in uV
                sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
                sb.hTextY_Pos= [-1000,-1000];
                
                fileName1 = fileName;
                if  contains(fileName,'_') %remove underscores for title
                    fileName1(strfind(fileName1,'_'))=' ';
                    fileName1=strcat('{',fileName1(1:end-4),'}');
                    title({[int2str(length(samples_index)/fs),'{ s of recording from:}'],[fileName1],...
                        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(sb.XLen),'{ s}','{ vertical }',...
                        int2str(sb.YLen),'{ uV}',')'],[' ']});
                else
                    title({[int2str(length(samples_index)/fs),'{ s of recording from:}'],[fileName1],...
                        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(sb.XLen),'{ s}','{ vertical }',...
                        int2str(sb.YLen),'{ uV}',')'],[' ']});
                end
                
                %make title font smaller
                ax = gca;
                ax.TitleFontSizeMultiplier = 0.9;
                
                %save as PNG
                fprintf(strcat('\n','\n',fileName(1:end-4),'\n','\n',' saving fig...', '\n','\n'))
                saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'s_trace_spikes_marked.png'));
                close all
                
            end
            
            %% find spikiest 100 ms
            close all
            
            time_win        =  time_wins(2); % in ms
            
            %check if already done
            if ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'ms_trace_spikes_marked.png'))
                
                fprintf(strcat('\n \n',fileName(1:end-4),'\n',...
                    '_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'ms_trace_spikes_marked.png \n \n already done! \n \n'))
                
            else %if not already done, carry out below...
                
                %if we don't have the spike train because above figure was already
                %complete, get it now...
                if  ~exist('spikeTrain')
                    spike_suffix_index = strfind(fileName,'Spikes')-3; %-3 to remove _aS in aSpikes for example
                    spike_suffix = fileName(spike_suffix_index+1:end-4);
                    raw_dat_fileName = strcat(fileName(1:spike_suffix_index),'.mat');
                    
                    if  ~~exist(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'))
                        load(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'));
                        finalData =     filteredMatrix(:,electrode_to_plot(ei));
                        clear filteredMatrix
                        spikeTrain =    spikeMatrix(:,electrode_to_plot(ei));
                        threshold = threshold;
                    else
                        load(raw_dat_fileName);
                        [spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot(ei)),method, multiplier,L,refPeriod_ms);
                    end
                    disp('acquired filtered data')
                end
                
                
                if ~exist('fs') %correction if fs isn't loaded properly (sampling freq)
                    fs = 25000;
                end
                
                sample_win      = time_win * (fs/1000); % because sampling freq, fs, is in seconds
                slidage         = 10 * fs/1000; % fs = 1 s; so this slides every 5 ms
                num_windows     = (length(spikeTrain)- sample_win) / slidage;
                num_done        = 0;
                %progressbar
                for i = 1:num_windows %-1 so that i can start from 0
                    spike_sums(i)   = sum(spikeTrain(1+slidage*num_done:slidage*num_done+sample_win));
                    num_done        = num_done + 1;
                    %progressbar(i/num_windows)
                end
                most_spikes_index = find(spike_sums == max(spike_sums))-1; % -1 is necessary
                %because if the first window had the most spikes, we need to start from
                %one; if the 2nd window had the most spikes, we need to start from
                %1+slidage*1
                
                if  length(most_spikes_index) > 1
                    most_spikes_index = most_spikes_index(2); %if it's a draw just take second one
                else
                end
                
                %get index of samples to plot
                samples_index = 1+slidage*most_spikes_index:slidage*most_spikes_index+sample_win;
                
                %plot
                figure
                x1 = (1:length(samples_index)) / (fs/1000); %time vec in ms
                try
                    y1 = finalData(samples_index);
                catch
                    samples_index = 1+slidage*most_spikes_index:length(finalData);
                    x1 = (1:length(samples_index)) / (fs/1000); %time vec in ms
                    y1 = finalData(samples_index);
                end
                plot(x1,y1,'Color','k');
                aesthetics
                box off
                %ylim([min(y1) max(y1)])
                ylim([-50 50])
                hold on
                %y2 = mean(y1)+6*std(y1);
                y2 = 50;
                sp_times = find(spikeTrain(samples_index) == 1);
                plot(x1(sp_times),y2,'v','LineWidth',0.5,'MarkerFaceColor','r','MarkerEdgeColor','r');
                
                removeAxis
                sb=scalebar;
                sb.Position=[0.5*(100/30),-40];
                sb.XLen = 2.5*(100/30); % in ms
                sb.YLen = 5; % in uV
                sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
                sb.hTextY_Pos= [-1000,-1000];
                
                fileName1 = fileName;
                if  contains(fileName,'_') %remove underscores for title
                    fileName1(strfind(fileName1,'_'))=' ';
                    fileName1=strcat('{',fileName1(1:end-4),'}');
                    title({[int2str(length(samples_index)/fs*1000),'{ ms of recording from:}'],[fileName1],...
                        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(floor(sb.XLen)),'{ ms}','{ vertical }',...
                        int2str(sb.YLen),'{ uV}',')'],[' ']});
                else
                    title({[int2str(length(samples_index)/fs*1000),'{ ms of recording from:}'],[fileName1],...
                        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(floor(sb.XLen)),'{ ms}','{ vertical }',...
                        int2str(sb.YLen),'{ uV}',')'],[' ']});
                end
                
                %make title font smaller
                ax = gca;
                ax.TitleFontSizeMultiplier = 0.9;
                
                %save as PNG
                fprintf(strcat('\n','\n',fileName(1:end-4),'\n','\n',' saving fig...', '\n','\n'))
                saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'ms_trace_spikes_marked.png'));
                close all
                
            end
            
            %% spikiest 30 ms
            close all
            
            time_win        = time_wins(3); % in ms
            
            %check if already done
            if ~~exist(strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'ms_trace_spikes_marked.png'))
                
                fprintf(strcat('\n \n',fileName(1:end-4),'\n',...
                    '_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'ms_trace_spikes_marked.png \n \n already done! \n \n'))
                
            else %if not already done, carry out below...
                
                %if we don't have the spike train because above figure was already
                %complete, get it now...
                if  ~exist('spikeTrain')
                    spike_suffix_index = strfind(fileName,'Spikes')-3; %-3 to remove _aS in aSpikes for example
                    spike_suffix = fileName(spike_suffix_index+1:end-4);
                    raw_dat_fileName = strcat(fileName(1:spike_suffix_index),'.mat');
                    
                    if  ~~exist(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'))
                        load(strcat(raw_dat_fileName(1:end-4),'_Filtd.mat'));
                        finalData =     filteredMatrix(:,electrode_to_plot(ei));
                        clear filteredMatrix
                        spikeTrain =    spikeMatrix(:,electrode_to_plot(ei));
                        threshold = threshold;
                    else
                        load(raw_dat_fileName);
                        [spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot(ei)),method, multiplier,L,refPeriod_ms);
                    end
                    disp('acquired filtered data')
                end
                
                
                if ~exist('fs') %correction if fs isn't loaded properly (sampling freq)
                    fs = 25000;
                end
                
                sample_win      = time_win * (fs/1000); % because sampling freq, fs, is in seconds
                slidage         = 5 * fs/1000; % fs = 1 s; so this slides every 5 ms
                num_windows     = (length(spikeTrain)- sample_win) / slidage;
                num_done        = 0;
                %progressbar
                for i = 1:num_windows %-1 so that i can start from 0
                    spike_sums(i)   = sum(spikeTrain(1+slidage*num_done:slidage*num_done+sample_win));
                    num_done        = num_done + 1;
                    %progressbar(i/num_windows)
                end
                most_spikes_index = find(spike_sums == max(spike_sums))-1; % -1 is necessary
                %because if the first window had the most spikes, we need to start from
                %one; if the 2nd window had the most spikes, we need to start from
                %1+slidage*1
                
                if  length(most_spikes_index) > 1
                    most_spikes_index = most_spikes_index(2); %if it's a draw just take second one
                else
                end
                
                %get index of samples to plot
                samples_index = 1+slidage*most_spikes_index:slidage*most_spikes_index+sample_win;
                
                %plot
                figure
                x1 = (1:length(samples_index)) / (fs/1000); %time vec in ms
                y1 = finalData(samples_index);
                plot(x1,y1,'Color','k');
                aesthetics
                box off
                %ylim([min(y1) max(y1)])
                ylim([-50 50])
                hold on
                %y2 = mean(y1)+6*std(y1);
                y2 = 50;
                sp_times = find(spikeTrain(samples_index) == 1);
                plot(x1(sp_times),y2,'v','LineWidth',0.5,'MarkerFaceColor','r','MarkerEdgeColor','r');
                
                removeAxis
                sb=scalebar;
                sb.Position=[0.5,-40];
                sb.XLen = 2.5; % in ms
                sb.YLen = 5; % in uV
                sb.hTextX_Pos= [-1000,-1000]; %-100 n both to make it disappear off screen
                sb.hTextY_Pos= [-1000,-1000];
                
                fileName1 = fileName;
                if  contains(fileName,'_') %remove underscores for title
                    fileName1(strfind(fileName1,'_'))=' ';
                    fileName1=strcat('{',fileName1(1:end-4),'}');
                    title({[int2str(length(samples_index)/fs*1000),'{ ms of recording from:}'],[fileName1],...
                        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(sb.XLen),'{ ms}','{ vertical }',...
                        int2str(sb.YLen),'{ uV}',')'],[' ']});
                else
                    title({[int2str(length(samples_index)/fs*1000),'{ ms of recording from:}'],[fileName1],...
                        ['{channel (MCS xy coordinate): }',num2str(channels(electrode_to_plot(ei)))],['{ (scalebars: horizontal }',num2str(sb.XLen),'{ ms}','{ vertical }',...
                        int2str(sb.YLen),'{ uV}',')'],[' ']});
                end
                
                %make title font smaller
                ax = gca;
                ax.TitleFontSizeMultiplier = 0.9;
                
                %save as PNG
                fprintf(strcat('\n','\n',fileName(1:end-4),'\n','\n',' saving fig...', '\n','\n'))
                saveas(gcf,strcat(fileName(1:end-4),'_channelXY_',int2str(channels(electrode_to_plot(ei))),...
                    '_',num2str(time_win),'ms_trace_spikes_marked.png'));
                close all
            end
        end
        
    end
    
end

end
