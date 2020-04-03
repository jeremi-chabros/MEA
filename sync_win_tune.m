% adjacency matrix sync window tuning

% currently only works for one file and does not save the file

%% get adjacency matrices
close all; clearvars -except refPeriod_ms

files = dir('*200127*FTD*slice2_mSpikes_4.mat');                   % where your .mat files are
files(2) = dir('*200127*FTD*slice3_mSpikes_4.mat');                   % where your .mat files are

files = files(~contains({files.name}, 'TTX'));
files = files(~contains({files.name}, 'ttx'));
files = files(~contains({files.name}, 'stim'));
files = files(~contains({files.name}, 'edited'));
files = files(~contains({files.name}, '2min'));
files = files(~contains({files.name}, 'adjM'));
files = files(~contains({files.name}, '191210'));


%%%%%%% set these parameters first:
method = 'tileCoef';
sync_win_s = [2 3 4 5 6]; %synchroncity window in seconds; e.g. 1 is +/- 1 s is considered synchronous (2DeltaT = 2 s)
rec_length_s = 360;
fs = 25000;
rec_length_samps = fs * rec_length_s;

%%%% downsampling:
num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
sync_windows = sync_win_s * ds_rate; %downsample the sync_win

progressbar
for file = 1:length(files)
    count = 1;
    load(files(file).name)
    for s = sync_windows
        sampling_fr = 25000;
        
        %get spike info
        try
            spikeMat=full(mSpikes);
            disp('loaded mSpikes')
        catch
            spikeMat=full(cSpikes);
            disp('loaded cSpikes')
        end
        
        %         if ~exist('channels')
        %             channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
        %             channels=channels.channels;
        %             disp('channels loaded from dat file')
        %         end
        
        spikeCounts=sum(spikeMat);
        %remove ref channel spikes:
        spikeCounts(find(channels == 15)) = 0;        %cell containing n. spikes for each channel
        active_chanIndex=find(spikeCounts>=10);
        ActiveSpikeCounts=spikeCounts(active_chanIndex);  %spikes of only active channels ('active'= >9)
        
        %get adjM
        fprintf(strcat('\n','\n','getting adjacency matrix...','\n','\n'))
        sync_win = sync_windows(count);
        adjM = getAdjM(spikeMat, method, num_samples,sync_win); %0 means no downsampling; 0.05 is sync window in s
        
        adjM1= adjM-eye(size(adjM));
        adjM1(find(isnan(adjM1)))=0;
        active_adjM=adjM1(active_chanIndex,active_chanIndex); %for extracting
        %part of full adjmat but now just compute active adjmat
        
        %mean correlation of active channels
        meanSTTC(file,count)= round(sum(sum(active_adjM))/(length(active_adjM)*(length(active_adjM)-1)),3);
        SEM_STTC(file,count)= round(std(sum(active_adjM))/(sqrt(length(active_adjM)*(length(active_adjM)-1))),3);
        count = count + 1;
    end
    figure
    x = 1:length(sync_windows);
    y = meanSTTC(file,:);
    plot(x,y,'Marker','s','LineStyle','none','MarkerFaceColor','red','MarkerSize',10,'MarkerEdgeColor','red')
    hold on
    errorbar(meanSTTC(file,:),SEM_STTC(file,:),'Color','red','LineWidth',2)
    aesthetics
    ylim([0 1])
    ylabel('Array mean spike time tiling coefficient')
    xticks(1:length(sync_windows))
    xticklabels(sync_win_s*1000)
    xlabel('Synchronicity window (ms)')
    if sum(sync_win_s*1000 > 999) >= 1
        xticklabels(round(sync_win_s,1));
        xlabel('Synchronicity window (s)')
    end
    set(gca,'TickDir','out')
    
    set(gca, 'FontSize', 14)
    
    fileName1=files(file).name;
    if  contains(files(file).name,'_') %remove underscores for title
        fileName1(strfind(fileName1,'_'))=' ';
        %fileName1=strcat('{',fileName1,'}');
        title({strcat(fileName1(1:end-4),' adjM parameter tuning'),' '});
    else
        title({strcat(files(file).name(1:end-4),' adjM parameter tuning'),' '});
    end
    
    %make title font smaller
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.7;
    
    fprintf(strcat('\n','\n',files(file).name(1:end-4),' saving adjM...', '\n','\n'))
    saveas(gcf,strcat(files(file).name(1:end-4),'_adjM_param_tune.png'));
    close all; clear adjM fileName fileName1
    
    progressbar(file/length(files))
    
end

%% load adjMs and save as PNG
% clear files adjM
% files = dir('*FTD*mSpikes_3_adjM_0.175.mat');  % where your .mat files are
%
%
%
%     fileName    =   files(i).name;
%     if ~exist(strcat(fileName(1:end-4),'_connectivity_matrix_',num2str(sync_win/ds_rate),'.png'))
%         load(fileName);
%         figure; imagesc(adjM);
%         aesthetics
%         ylabel('Electrode')
%         xlabel('Electrode')
%         cb = colorbar;
%         % ylabel(cb, 'Spike count')
%         ylabel(cb, 'Correlation')
%         cb.TickDirection = 'out';
%         % cb.Ticks = 0:5; % for slice 5 specifically
%         set(gca,'TickDir','out');
%         %cb.Location = 'Southoutside';
%         cb.Box = 'off';
%         set(gca, 'FontSize', 14)
%         caxis([0,1])
%         %determine file name and title of fig
%         fileName1=files(i).name;
%         if  contains(files(i).name,'_') %remove underscores for title
%             fileName1(strfind(fileName1,'_'))=' ';
%             %fileName1=strcat('{',fileName1,'}');
%             title({strcat(fileName1(1:end-4),' connectivity matrix'),' '});
%         else
%             title({strcat(files(i).name(1:end-4),' connectivity matrix'),' '});
%         end
%         %make title font smaller
%         ax = gca;
%         ax.TitleFontSizeMultiplier = 0.7;
%         %save as PNG
%         fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving adjM...', '\n','\n'))
%         saveas(gcf,strcat(fileName(1:end-4),'_connectivity_matrix_',num2str(sync_win),'.png'));
%         close all; clear adjM fileName fileName1
%     else
%     end
% end
%
%     end
