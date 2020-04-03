
%% PV-ArchT 
% figure;histogram(output(2).spikes)
% 
% figure;histogram(output(8).spikes)
load('cwt0_PVArchT.mat')

h = figure;histogram(output(8).spikes,'BinWidth',50);
hold on;histogram(output(2).spikes,'BinWidth',50,'BinLimits',[500 4000])
hold on;histogram(output(5).spikes,'BinWidth',50)
xlim([0 3000])
% was going to change y lim to 59
% ed_y = yticks; 
% ed_y(end) = 59
aesthetics
ylabel('Frequency')
xlabel('Spike count')
legend('Baseline','Light on','Light off')
set(gca,'fontname','arial','fontsize',22)
title({'PV-ArchT (200219 2C DIV17 cSpikes L0)',' '},'fontsize',12)
yticks([0 15 30 45 59])
yticklabels([0 15 30 45 60])
set(gcf,'Position',[680 593 425 385])

saveas(h,'PAT200219_2C_DIV17_cSpikes_L0.mat_hist_spikecounts.png')
close(h)

%% PV-Ai32
load('cwt0_PVAi32.mat')
% 190904 1E DIV26 cSpikes L0
h = figure; h1 = histogram(log10(output(2).spikes));
hold on;histogram(log10(output(1).spikes))
aesthetics
ylabel('Frequency')
xlabel('Spike count')
legend('Baseline','Light on')
set(gca,'fontname','arial','fontsize',22)
title({'PV-Ai32 (190904 1E DIV26 cSpikes L0)',' '},'fontsize',12)
set(gcf,'Position',[680 593 425 385])
xticklabels(10.^xticks)
yts = yticks;
xts = xlim;
saveas(h,'190904_1E_DIV26_cSpikes_L0.mat_hist_spikecounts.png')
close(h)

% 190904 2E DIV26 cSpikes L0
h = figure;h1 = histogram(log10(output(4).spikes),'BinWidth',0.5);
hold on;histogram(log10(output(3).spikes),'BinWidth',0.5)
aesthetics
ylabel('Frequency')
xlabel('Spike count')
legend('Baseline','Light on')
set(gca,'fontname','arial','fontsize',22)
title({'PV-Ai32 (190904 2E DIV26 cSpikes L0)',' '},'fontsize',12)
set(gcf,'Position',[680 593 425 385])
ylim([min(yts) max(yts)])
xlim(xts)
xticklabels(10.^xticks)
saveas(h,'190904_2E_DIV26_cSpikes_L0.mat_hist_spikecounts.png')
close(h)