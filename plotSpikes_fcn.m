%% plot filtered trace of spikes overlaid and peaks aligned 
function plotSpikes_fcn(channel,file,method,parameter)
clear all
spikeTrain=zeros(1500000,1);
finalData=zeros(1500000,1);
dat=zeros(1500000,60);
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
load('200114_FTDOrg_GrpA_5A_Slice6.mat')
electrode_to_plot=53; %not MCS ID but which column the desired E is in
multiplier=5;
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'Manuel', multiplier,0);
%dat=dat(:,3); 
%[spikeTrain, finalData, threshold] = detectSpikes(dat, 'Manuel', 4, 0);
%spikeTrain=load('MPT190515_4B_DIV21_cSpikes_L0.mat');
%spikeTrain=spikeTrain.cSpikes;
%spikeTrain=full(spikeTrain);
%spikeTrain=spikeTrain(:,27);

sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=50;
%added correction if num spikes is fewer than desired number:
if sum(spikeTrain) < n_spikes_to_plot
    n_spikes_to_plot = sum(spikeTrain);
end
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
ylim([-90 50])
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

%find(v_trace_ch_47(sp_times(203):sp_times(203)+25)==min(v_trace_ch_47(sp_times(203):sp_times(203)+25));

%% wavelet
%load
clear all
%spikeTrain=zeros(18000000,1);
%finalData=zeros(18000000,1);
%dat=zeros(18000000,60);
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
filename='SMPT190923_2B_DIV21.mat';
load(filename)
%% run cwt
%dat=dat(:,27); 
electrode_to_plot=15; %not MCS ID but which column the desired E is in
L=0.1254;
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'cwt', 0,L);
%[spikeTrain, finalData, threshold] = detectSpikes(dat, 'cwt', 0,L);
%spikeTrain=load('MPT190515_4B_DIV21_cSpikes_L0.mat');
%spikeTrain=spikeTrain.cSpikes;
%spikeTrain=full(spikeTrain);
%spikeTrain=spikeTrain(:,27);
%% plot cwt
sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=1;
%n_spikes_to_plot=length(sp_times);
for i=1:n_spikes_to_plot
    sp_peak_time=find(finalData(sp_times(i):sp_times(i)+25)==min(finalData(sp_times(i):sp_times(i)+25)));%note changed min to max as with filtered data spike is +ve
    %plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25));
    plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25),...
        'Color',[0.5 0.5 0.5],'LineWidth',0.1); %all grey
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
ylim([-200 100])
yticks(linspace(-200,100,(300/50)+1)); %last # is num out (should be 'num from' minus 'num to' plus one)
ylabel('filtered signal (\muV)');
set(gca,'fontsize',16)
%threshvec=ones(length(ave_trace))*threshold;
%hold on
%plot(1:length(ave_trace),threshvec,'LineStyle','--','Color','b','LineWidth',3)
%hold off
%calculate the average and overlay that onto the figure as a red line
%surrounded by white
title({[int2str(n_spikes_to_plot),'{ template spikes (L = }',num2str(L),')'],['{ from }',filename],...
    [' channel: ',int2str(channels(electrode_to_plot))]})
clear ave_trace

%% plot abs method
clear all
spikeTrain=zeros(1500000,1);
finalData=zeros(1500000,1);
dat=zeros(1500000,60);
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
load('200114_FTDOrg_GrpA_5A_Slice1.mat')
electrode_to_plot=40; %not MCS ID but which column the desired E is in
multiplier=-27;
[spikeTrain, finalData, threshold] = detectSpikes(dat(:,electrode_to_plot), 'abs', multiplier,0);
%dat=dat(:,3); 
%[spikeTrain, finalData, threshold] = detectSpikes(dat, 'Manuel', 4, 0);
%spikeTrain=load('MPT190515_4B_DIV21_cSpikes_L0.mat');
%spikeTrain=spikeTrain.cSpikes;
%spikeTrain=full(spikeTrain);
%spikeTrain=spikeTrain(:,27);

sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=40;
%added correction if num spikes is fewer than desired number:
if sum(spikeTrain) < n_spikes_to_plot
    n_spikes_to_plot = sum(spikeTrain);
end
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
ylim([-90 50])
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

%find(v_trace_ch_47(sp_times(203):sp_times(203)+25)==min(v_trace_ch_47(sp_times(203):sp_times(203)+25));

end