%% plot filtered trace of spikes overlaid and peaks aligned 
clear all
spikeTrain=zeros(18000000,1);
finalData=zeros(18000000,1);
dat=zeros(18000000,60);
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
load('MPT190515_4B_DIV21.mat')
dat=dat(:,27); 
[spikeTrain, finalData, threshold] = detectSpikes(dat, 'Manuel', 5, 0);
%spikeTrain=load('MPT190515_4B_DIV21_cSpikes_L0.mat');
%spikeTrain=spikeTrain.cSpikes;
%spikeTrain=full(spikeTrain);
%spikeTrain=spikeTrain(:,27);

sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=50;
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

%change axes to voltage normalised
%change x axis to time rather than samples
xticks(linspace(0,50,3));
xticklabels(linspace(0,50,3)/25);
xlabel('time (ms)');
ylim([-80 80])
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
clear all
spikeTrain=zeros(18000000,1);
finalData=zeros(18000000,1);
dat=zeros(18000000,60);
%cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.MAT_Files_Whole'
load('MPT190515_4B_DIV21.mat')
dat=dat(:,27); 
L=0;
[spikeTrain, finalData, threshold] = detectSpikes(dat, 'cwt', 0,L);
%spikeTrain=load('MPT190515_4B_DIV21_cSpikes_L0.mat');
%spikeTrain=spikeTrain.cSpikes;
%spikeTrain=full(spikeTrain);
%spikeTrain=spikeTrain(:,27);

sp_times=find(spikeTrain==1);
figure
n_spikes_to_plot=50;
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


%change axes to voltage normalised
%change x axis to time rather than samples
xticks(linspace(0,50,3));
xticklabels(linspace(0,50,3)/25);
xlabel('time (ms)');
ylim([-80 80])
yticks(linspace(-80,80,(160/20)+1));
ylabel('filtered signal (\muV)');
set(gca,'fontsize',16)
%threshvec=ones(length(ave_trace))*threshold;
%hold on
%plot(1:length(ave_trace),threshvec,'LineStyle','--','Color','b','LineWidth',3)
%hold off
%calculate the average and overlay that onto the figure as a red line
%surrounded by white
clear ave_trace