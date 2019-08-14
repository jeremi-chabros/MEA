% plot spike rate over development

%open directory with spike matrices
%cd 'D:\MECP2_2019_AD\Data_To_Use\4.0.SpikeMatrices'
    files = dir('*cSpikes*');  % where your .mat files are 
    %go through each file and sum spikes depending on DIV
    
    %initialise
    DIV07Spikes=[];
    DIV14Spikes=[];
    DIV21Spikes=[];
    DIV28Spikes=[];
    DIV35Spikes=[];
    
for i=1:length(files)
    load(files(i).name)
    Spikes = full(cSpikes);
    if files(i).name(14:18)=='DIV07'
        DIV07Spikes(i)=full(sum(sum(Spikes)));
    elseif files(i).name(14:18)=='DIV14'
        DIV14Spikes(i)=full(sum(sum(Spikes)));
    elseif files(i).name(14:18)=='DIV21'
        DIV21Spikes(i)=full(sum(sum(Spikes)));
    elseif files(i).name(14:18)=='DIV28'
        DIV28Spikes(i)=full(sum(sum(Spikes)));
    elseif files(i).name(14:18)=='DIV35'
        DIV35Spikes(i)=full(sum(sum(Spikes)));    
    else
    end             
end 

    DIV07Spikes=((sum(DIV07Spikes))/sum(DIV07Spikes>0))/60;
    DIV14Spikes=((sum(DIV14Spikes))/sum(DIV14Spikes>0))/60;
    DIV21Spikes=((sum(DIV21Spikes))/sum(DIV21Spikes>0))/60;
    DIV28Spikes=((sum(DIV28Spikes))/sum(DIV28Spikes>0))/60;
    DIV35Spikes=((sum(DIV35Spikes))/sum(DIV35Spikes>0))/60;
    
    S_over_dev=round([DIV07Spikes,DIV14Spikes,DIV21Spikes,DIV28Spikes,DIV35Spikes]);
    
%plot
figure
bar(S_over_dev);
xticklabels([07,14,21,28,35])
xlabel('Days in vitro')
ylabel('average number of spikes across all electrodes of culture')
box off
title('how spike count changes over development')
aesthetics
%culture
% culture='MPT010119_4B';

%load file
%load('MPT070119_4B_DIV7_info.mat');


%% plot distribution of spike counts on DIV14
DIV14Spikes_dis=zeros(1,60);
for i=1:length(files)
    load(files(i).name)
    disp(files(i).name)
    if files(i).name(14:18)=='DIV14'
        DIV14Spikes_d=full(sum(Spikes));
        DIV14Spikes_dis(i,:)=DIV14Spikes_d;
    else
    end             
end 
%take only non zero elements
DIV14Spikes_dis=DIV14Spikes_dis(DIV14Spikes_dis>0);

figure 
hist(DIV14Spikes_dis)
title('distribution of spike counts')
box off

figure 
hist(log(DIV14Spikes_dis))
title('distribution of ln spike counts (should be 2 peaks)')
box off

%% plot spikes overlaid and peaks aligned 

%cd 'D:\MECP2_2019_AD\Data_To_Use\4.0.SpikeMatrices'
%load('MPT190403_6C_DIV14_info.mat')
sp_train_ch_47=full(tSpikes(:,1));
cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
load('MPT190403_6C_DIV14.mat')
v_trace_ch_47=dat(:,1);

sp_times=find(sp_train_ch_47==1);
figure
for i=1:length(sp_times)
    sp_peak_time=find(v_trace_ch_47(sp_times(i):sp_times(i)+25)==min(v_trace_ch_47(sp_times(i):sp_times(i)+25)));
    plot(v_trace_ch_47((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25))
    hold on
        ave_trace(:,i)=v_trace_ch_47((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25);

end
ave_trace=(sum(ave_trace'))/(length(ave_trace(1,:)));
plot(ave_trace,'LineWidth',8,'Color','w');
plot(ave_trace,'LineWidth',3,'Color','r');
hold off
%aesthetics
box off
clear ave_trace
%change axes to voltage normalised
%change x axis to time rather than samples
xticks(linspace(0,60,7));
xticklabels(linspace(0,60,7)/25);
xlabel('time (ms)');
ylabel('uV');
%calculate the average and overlay that onto the figure as a red line
%surrounded by white


%find(v_trace_ch_47(sp_times(203):sp_times(203)+25)==min(v_trace_ch_47(sp_times(203):sp_times(203)+25));


%% plot filtered trace of spikes overlaid and peaks aligned 
clear all
spikeTrain=zeros(18000000,1);
finalData=zeros(18000000,1);
dat=zeros(18000000,60);
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
load('MPT190515_4B_DIV21.mat')
dat=dat(:,27); 
[spikeTrain, finalData, threshold] = detectSpikes(dat, 'Manuel', 5, 0);
spikeTrain=load('MPT190515_4B_DIV21_cSpikes_L0.mat');
spikeTrain=spikeTrain.cSpikes;
spikeTrain=full(spikeTrain);
spikeTrain=spikeTrain(:,27);

sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=5;
%n_spikes_to_plot=length(sp_times);
for i=1:n_spikes_to_plot
    sp_peak_time=find(finalData(sp_times(i):sp_times(i)+25)==min(finalData(sp_times(i):sp_times(i)+25)));%note changed min to max as with filtered data spike is +ve
    plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25));
    hold on
    ave_trace(:,i)=finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25);
    
end
ave_trace=(sum(ave_trace'))/(length(ave_trace(1,:)));
plot(ave_trace,'LineWidth',8,'Color','w');
plot(ave_trace,'LineWidth',3,'Color','r');
hold off
aesthetics
box off
clear ave_trace
%change axes to voltage normalised
%change x axis to time rather than samples
xticks(linspace(0,60,7));
xticklabels(linspace(0,60,7)/25);
xlabel('time (ms)');
ylabel('filtered signal (\muV)');
%calculate the average and overlay that onto the figure as a red line
%surrounded by white


%find(v_trace_ch_47(sp_times(203):sp_times(203)+25)==min(v_trace_ch_47(sp_times(203):sp_times(203)+25));

%% same again but use downsample trace?

