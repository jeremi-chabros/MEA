% correct adjM remove edges from ref channel

clear all
filename='SMPT190923_2B_DIV21_cSpikes_L0.1254_adjM_reord.mat';
channelID=15; %channel that you want to ground (use MCS coordinate)
load(filename);
ref=find(new_channels==15);
%for i =1:size(adjM,1)
%grd_channels(i,1:size(adjM,1))=adjM(i,:)==adjM(ref,:);
%end
%channels not identical to ref


%remove 1s then add IdMat back in
adjM(find(adjM==1))=0;
adjM=adjM+flip(eye(size(adjM)));
figure;imagesc(adjM);
yticks([0 10 20 30 40 50 60])
yticklabels(flip(yticklabels))
yticks([1 10 20 30 40 50 60])


fileName = strcat(filename(1:end-4), '_2.mat'); 
        % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
        save(fileName, 'adjM','channels','new_channels');