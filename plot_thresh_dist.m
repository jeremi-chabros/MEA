%code used to generate distribution of thesholds

load('190814_FTD_slice4_mSpikes_3.mat')
channel(49)
channels(49)
full(sum(mSpikes))
hist(thresholds)
aesthetics
ylab('Frequency')
ylabel('Frequency')
xlabel('Threshold ('uV')')
xlabel('Threshold (uV)')
mean(thresholds)
m=mean(thresholds)
figure
hist(thresholds)
xline(m)
xline(m,'color','red')
xline(m,'color','red','linewidth','2')
xline(m,'color','red','linewidth',2)