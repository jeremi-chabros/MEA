
% script for comparing groups

%get data
cd 'C:\Users\alexd\Dropbox (Cambridge University)\NOG MEA Data\MEA2100 Organoid\Analysis\Unsorted'
load('manuel_onefile_3_organoid_stats_ages.mat')
%back to scripts directory
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'

%% betweenness centrality
clear groups
for i = 1:length(output)
    groups(i) = output(i).grp;
end
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).BC_n;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).BC_n;
end

close all; notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('normalised betweenness centraility')

%% degree
figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).meanDegree;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).meanDegree;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('average node degree')


%% nrmlsd_cnnctvty

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).nrmlsd_cnnctvty;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).nrmlsd_cnnctvty;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('normalised average functional connectivity')

%% meanSTTC

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).meanSTTC;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).meanSTTC;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('average functional connectivity')

%% firing rate


figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).mean_FR;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).mean_FR;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('average firing rate (spikes/s)')

%% netw_density

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).netw_density;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).netw_density;
end

notBoxPlot([dataE*100],5,'style','sdline')
hold on
notBoxPlot([dataD*100],4,'style','sdline')
notBoxPlot([dataC*100],3,'style','sdline')
notBoxPlot([dataB*100],2,'style','sdline')
notBoxPlot([dataA*100],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Network density (%)')

%% netw_size

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).netw_size;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).netw_size;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Network size (# nodes)')

%% CC


figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).CC;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).CC;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Normalised clustering')

%% PL

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).PL;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).PL;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Normalised characteristic path length')
%% SW

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).SW;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).SW;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).SW;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).SW;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).SW;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Small world coefficient')

%% RC

figure
ind = find(groups =='E');

for i = 1:length(ind)
    dataE(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='D');
for i = 1:length(ind)
    dataD(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='C');
for i = 1:length(ind)
    dataC(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='B');
for i = 1:length(ind)
    dataB(i) = output(ind(i)).RC;
end

clear ind; ind = find(groups =='A');
for i = 1:length(ind)
    dataA(i) = output(ind(i)).RC;
end

notBoxPlot([dataE],5,'style','sdline')
hold on
notBoxPlot([dataD],4,'style','sdline')
notBoxPlot([dataC],3,'style','sdline')
notBoxPlot([dataB],2,'style','sdline')
notBoxPlot([dataA],1,'style','sdline')

aesthetics
xticks(1:5)
xticklabels({'CS30','EpiC','WTsli42','D','CS30 older'})
ylabel('Rich club coefficient')

%%  +AGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% firing rate ~ group * age

clear groups
for i = 1:length(output)
    
    groups{i} = strcat(num2str(output(i).age_DIV),output(i).grp);
    
end

group_names = unique(groups);

clear data
for j = 1:length(group_names)
    
    clear ind
    ind = strcmp(groups,group_names{j});
    ind = find(ind == 1);
    
    for i = 1:length(ind)
        data_temp(i,1) = output(ind(i)).mean_FR;
    end
    
    data{j} = data_temp;
    clear data_temp

end

% figure
% cols = [linspace(0,1,length(group_names))',linspace(1,0,length(group_names))',linspace(0,1,length(group_names))'];
% perm = [6     3     7     1     5    10     9    12     8     2    11     4];
% cols(:,1:3) = cols(perm,1:3);
% 
% rand
% for i = 1:length(group_names)
%     
%     grp = group_names{i};
%     b = notBoxPlot(data{i},str2num(grp(1:3)),'style','sdline');
%     b.data.MarkerFaceColor = cols(i,1:3);
%     hold on
%     clear grp
%     
% end
% aesthetics
% xlim([140 250])


figure
grp = group_names{1};
b = notBoxPlot(data{1},str2num(grp(1:3)),'style','sdline');
b.mu.Color = [1 1 1];
b.semPtch.FaceColor = [1 1 1]
b.sd.Color = [1 1 1];
b.semPtch.EdgeColor = [1 1 1];

hold on
b = notBoxPlot(data{2},str2num(grp(1:3))-0.2,'style','sdline');
b.data.Marker = 'd';
b.mu.Color = [1 1 1];
b.semPtch.FaceColor = [1 1 1]
b.sd.Color = [1 1 1];
b.semPtch.EdgeColor = [1 1 1];

hold on
b = notBoxPlot(data{3},str2num(grp(1:3))+0.2,'style','sdline');
b.data.Marker = 's';
b.mu.Color = [1 1 1];
b.semPtch.FaceColor = [1 1 1]
b.sd.Color = [1 1 1];
b.semPtch.EdgeColor = [1 1 1];















