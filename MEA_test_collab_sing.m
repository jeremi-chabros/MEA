% analysis for MEA testing collaboration with Andras, Sagnik, Ilaria
clear all
%% get files and save spikes
clear all; close all
% inputs:
data_and_scripts_dir = 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output';
files = dir('*edited*.mat*');                   % where your .mat files are 
files = files(~contains({files.name}, 'Spikes'));
%files = files(~contains({files.name}, 'slice1'));
%files = files(~contains({files.name}, '_2'));
meths   ={'Manuel'};
%meths   ={'cwt'};                              % one method only for speed
%params  =[4,5,6; ...                            % Threshold  
%    -0.1254,0,0.1254]';                         % L parameter
params  = [6, 0];                       %one threshold, one template
% code:
disp(strcat({'number of files to do = '},num2str(length(files))))
progressbar('files','parameters','methods')
for m=1:length(meths)
    method              =       meths{m}
    
    for p = 1:size(params,1);
        L               =       params(p+size(params,1))
        multiplier      =       params(p)
        
        batchGetSpike_function(data_and_scripts_dir,files,method,multiplier,L)
        progressbar([],p/size(params,1),[]);       %update parameter progress
    end
    progressbar([],[],m/length(meths));            %update % methods done
end
    
%% creater rasters and save
close all; clear all
%inputs:
files = dir('*191210*slice1*g04*mSpikes_6*.mat*');                   % where your .mat files are 
%files = dir('*FTD*Group*Spikes*.mat*');                   % where your .mat files are 
%files = files(~contains({files.name}, 'slice1'));
%files = files(~contains({files.name}, 'Spikes'));
%files = files(~contains({files.name}, '_4'));       % remove ones already done
files = files(~contains({files.name}, 'adjM'));
files = files(~contains({files.name}, '8_m'));

%code:
for i=1:length(files)
    fileName    =   files(i).name
    load(fileName)
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
    if  floor(length(spikeMatrix)/fs)~=length(spikeMatrix)/fs; 
        %calculate number of samples to subtract to make 
        %length of rec in s a whole number
        n2del = fs*(length(spikeMatrix)/fs - floor(length(spikeMatrix)/fs));
        spikeMatrix=spikeMatrix(1:length(spikeMatrix)-(n2del-1),:);
    else
    end

recordDuration = round(length(spikeMatrix)); %in samples
%downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
downFactor = 25000;
downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration/downFactor); 
new_fs = fs/downFactor;

figure
imagesc(downSpikeMatrix')
%imagesc(log(downSpikeMatrix'))
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
%caxis([0,3])
yticks([1, 10:10:60])
%cb.Limits = [0 4]; %only changes max number on scale bar doesnt change
%scale bar color
caxis([0,4]); 

        fileName1=files(i).name;
    if  contains(files(i).name,'_') %remove underscores for title
        fileName1(strfind(fileName1,'_'))=' ';
        %fileName1=strcat('{',fileName1,'}');
        title(strcat(fileName1(1:end-4),' Raster'));
    else
        title(strcat(files(i).name(1:end-4),' Raster'));
    end

%save raster as PNG
    fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving raster...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_Raster.png'));    
    close all; clear cSpikes mSpikes spikeMatrix fs downSpikeMatrix fileName fileName1
end

%% plot before and after ttx

%beforeLabel     =    'before TTX';   
%afterLabel      =    'after TTX';
%files = dir('*FTD*Group*Spikes*.mat*');                  
%files = files(~contains({files.name}, 'slice1'));


%beforeFile      =   
%afterFile       =   
%bar_chart_organoid_basic_fcn(beforeFile,afterFile,beforeLabel,afterLabel)

%% get adjacency matrices and save as PNG
clear all
method = 'tileCoef';
files = dir('*191210*slice1*Spikes_5*.mat*');                   % where your .mat files are 
files = files(~contains({files.name}, 'adjM'));

batch_getAdj_fcn(method,files);

% load adjMs and save as PNG
clear files adjM
files = dir('*191210*slice1*Spikes_5*adjM*.mat*');  % where your .mat files are

for i=1:length(files)
    fileName    =   files(i).name;
 if ~exist(strcat(fileName(1:end-4),'_connectivity_matrix.png'))
    load(fileName);
    figure; imagesc(adjM);
    aesthetics
    ylabel('Electrode') 
    xlabel('Electrode')
    cb = colorbar;
    % ylabel(cb, 'Spike count')
    ylabel(cb, 'Correlation') 
    cb.TickDirection = 'out';
    % cb.Ticks = 0:5; % for slice 5 specifically
    set(gca,'TickDir','out'); 
    %cb.Location = 'Southoutside';
    cb.Box = 'off';
    set(gca, 'FontSize', 14)
    
    %determine file name and title of fig
    fileName1=files(i).name;
    if  contains(files(i).name,'_') %remove underscores for title
        fileName1(strfind(fileName1,'_'))=' ';
        %fileName1=strcat('{',fileName1,'}');
        title(strcat(fileName1(1:end-4),' connectivity_matrix'));
    else
        title(strcat(files(i).name(1:end-4),' connectivity_matrix'));
    end
    
    %save as PNG
    fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving adjM...', '\n','\n'))
    saveas(gcf,strcat(fileName(1:end-4),'_connectivity_matrix.png'));    
    close all; clear adjM fileName fileName1
 else
 end    
end


%% reorder adjMs in order to use network plot in R
clear files
files = dir('*191210*slice1*Spikes_5*ADJm*.mat*');  % where your .mat files are
files = files(~contains({files.name}, 'reord'));
reorder_adjM_fcn(files);

%% to add part that calls R script, plots network then save as PNG

%% save heatMaps
clear files
files = dir('191210_FTD_slice1_DIV_g04_2018*_mSpikes_6.mat*');  % where your .mat files are
%files = files(~contains({files.name}, 'adjM'));
option = 'logc'; %logc, count or rate (gets capped around 10,000 spikes)
batch_getHeatMaps_fcn(files,option)

%% spikes overlay
% identify channels with most spikes
% overlay 50 spikes (or if less than 50, overlay these and note no. spikes
% in the title
% could create a function that overlays spikes in every channel (max of 50
% spikes)