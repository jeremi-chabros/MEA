
beforeFile='190830_slice1stim5_mSpikes_4_org.mat';
load(beforeFile)
try
    pre_spikemat=full(cSpikes);
catch
    pre_spikemat=full(mSpikes);
end

afterFile='190830_slice1_mSpikes_4_org.mat';
load(afterFile)
try
    pst_spikemat=full(cSpikes);
catch
    pst_spikemat=full(mSpikes);
end

load('190830_slice2_prestim_cSpikes_L-0.1254.mat')
pre_spikemat=full(cSpikes);
load('190830_slice2stim10002_cSpikes_L-0.1254.mat')
pst_spikemat=full(cSpikes);

%if pre is longer than post...
rec_size_ratio=length(pre_spikemat)/length(pst_spikemat);

pre=round(sum(pre_spikemat)/rec_size_ratio);
a_pre=pre(find(pre>=1));

pst=sum(pst_spikemat);
pst(15)=0; %remove ref channel noise  (there were 14 spikes)
a_pst=pst(find(pst>=1));

pre_mean1=mean(a_pre);
pre_stderror1= std(a_pre) / sqrt( length(a_pre))
pre_medi1=median(a_pre);
pre_acti1=length(a_pre);

pst_mean1=mean(a_pst);
pst_stderror1= std(a_pst) / sqrt( length(a_pre))
pst_medi1=median(a_pst);
pst_acti1=length(a_pst);

y=[pre_mean1 pst_mean1];
figure
bar(y)
aesthetics
xticklabels({'No stim.'; 'Stim'})
ylabel('Mean spike count')
title('30/08/19 slice two, stim protocol 3')
ylim([0 100])

x = 1:2;
data = y';
%errhigh = [y+[pre_stderror1 pst_stderror1]]';
%errlow  = [y-[pre_stderror1 pst_stderror1]]';
errhigh = [pre_stderror1 pst_stderror1]';
errlow  = [pre_stderror1 pst_stderror1]';

hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

figure
y=[pre_acti1 pst_acti1];
bar(y)
aesthetics
xticklabels({'No stim.'; 'Stim'})
ylabel('number of channels with spikes')
title('30/08/19 slice two, stim protocol 3')
ylim([0 60])

x = 1:2;
data = y';
%errhigh = [y+[pre_stderror1 pst_stderror1]]';
%errlow  = [y-[pre_stderror1 pst_stderror1]]';
errhigh = [pre_stderror1 pst_stderror1]';
errlow  = [pre_stderror1 pst_stderror1]';

hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off




