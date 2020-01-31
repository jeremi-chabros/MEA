function bar_chart_organoid_basic_fcn(beforeFile,afterFile)

%beforeFile='191209_FTD_slice1_GroupB_1_mSpikes_5.mat';
load(beforeFile)
try
    pre_spikemat=full(cSpikes);
catch
    pre_spikemat=full(mSpikes);
end

afterFile='191209_FTD_slice1_GroupB_2_ttx_mSpikes_5.mat';
load(afterFile)
try
    pst_spikemat=full(cSpikes);
catch
    pst_spikemat=full(mSpikes);
end
%remove spikes from ref
pre_spikemat(:,find(channels==15))=zeros(size(pre_spikemat(:,find(channels==15))));
pst_spikemat(:,find(channels==15))=zeros(size(pst_spikemat(:,find(channels==15))));

%if pre is longer than post...
rec_size_ratio=length(pre_spikemat)/length(pst_spikemat);

pre=round(sum(pre_spikemat)/rec_size_ratio);
a_pre=pre(find(pre>=10));

pst=sum(pst_spikemat);
pst(15)=0; %remove ref channel noise  (there were 14 spikes)
a_pst=pst(find(pst>=10));

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
%xticklabels({'No stim.'; 'Stim'})
xticklabels({beforeLabel; afterLabel})
ylabel('Mean spike count')
ylim([0 200])
if strfind(afterFile,'_') %remove underscores for title
    afterFile1=afterFile;
    afterFile1(strfind(afterFile,'_'))=' ';
    afterFile1=strcat('{',afterFile1,'}');
    title(afterFile1);
else
    title(afterFile);
end

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
%xticklabels({'No stim.'; 'Stim'})
xticklabels({beforeLabel; afterLabel})
ylabel('number of channels with >10 spikes')
ylim([0 60])
if strfind(afterFile,'_') %remove underscores for title
    afterFile1=afterFile;
    afterFile1(strfind(afterFile,'_'))=' ';
    afterFile1=strcat('{',afterFile1,'}');
    title(afterFile1);
else
    title(afterFile);
end

x = 1:2;
data = y'; 
%don't need error bars for n active channels as not an average
%errhigh = [y+[pre_stderror1 pst_stderror1]]';
%errlow  = [y-[pre_stderror1 pst_stderror1]]';
%errhigh = [pre_stderror1 pst_stderror1]';
%errlow  = [pre_stderror1 pst_stderror1]';

hold on

%er = errorbar(x,data,errlow,errhigh);    
%er.Color = [0 0 0];                            
%er.LineStyle = 'none';  

hold off

end



