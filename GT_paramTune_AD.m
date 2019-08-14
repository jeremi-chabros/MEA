% tune graph parameters

%load data
spikeMat=full(cSpikes);


%get adjM
method = 'tileCoef';
adjM = getAdjM(spikeMat, method, 0,0.05); %0 means no downsampling; 0.05 is sync window in s
adjM1= adjM-eye(size(adjM));
adjM1(find(isnan(adjM1)))=0;


 thresholds = [prctile(adjM1(:), 90), prctile(adjM1(:), 92.5),...
     prctile(adjM1(:), 95), prctile(adjM1(:), 97.5)]; % get the 95th percentile 

 %% all nodes (not just active) 


 for th = 1:length(thresholds)
     edges=adjM1;
edges(edges < thresholds(th)) = 0.0001; 
buEdges=edges;
buEdges(find(buEdges==0.0001))=0;
buEdges(find(buEdges>0.0001))=1;
density(th)=sum(sum(buEdges))/(length(buEdges)*(length(buEdges)-1));

thresh_P = [90 92.5 95 97.5];

%% active only
%active_chanIndex=find(spikeCounts>9);
%active_adjM=adjM1(active_chanIndex,active_chanIndex);
%active_edges=buEdges(active_chanIndex,active_chanIndex);
%density(th)=sum(sum(active_edges))/(length(active_edges)*(length(active_edges)-1));
%degree(th)=mean(sum(active_edges));

 end
 
  figure
 subplot(1,2,1)
     plot(thresh_P, density);
    xlabel('Threshold as a percentile') 
    ylabel('network density')
    lineThickness(1.5) 
    aesthetics
%% active
active_chanIndex=find(spikeCounts>9);
active_adjM=adjM1(active_chanIndex,active_chanIndex);

 %thresholds = [prctile(active_adjM(:), 10),prctile(active_adjM(:), 20),...
  %   prctile(active_adjM(:), 30),prctile(active_adjM(:), 40),...
   %  prctile(active_adjM(:), 50), prctile(active_adjM(:), 60),...
    % prctile(active_adjM(:), 70), prctile(active_adjM(:), 80),...
     %prctile(active_adjM(:), 90),prctile(active_adjM(:), 100)]; % get the 95th percentile 
thresholds = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

 for th = 1:length(thresholds)
     edges=active_adjM;
edges(edges < thresholds(th)) = 0.0001; 
buEdges=edges;
buEdges(find(buEdges==0.0001))=0;
buEdges(find(buEdges>0.0001))=1;
density(th)=sum(sum(buEdges))/(length(buEdges)*(length(buEdges)-1));

thresh_P = [90 92.5 95 97.5];

%% active only

%density(th)=sum(sum(active_edges))/(length(active_edges)*(length(active_edges)-1));
%degree(th)=mean(sum(active_edges));

 end
 
  figure
 subplot(1,2,1)
     scatter(thresholds, density,'+','r');
    xlabel('STTC threshold') 
    ylabel('network density')
    lineThickness(1.5) 
    aesthetics
    ylim([0 1])
 
 %% degree
 thresholds = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

 for th = 1:length(thresholds)
     edges=active_adjM;
edges(edges < thresholds(th)) = 0.0001; 
buEdges=edges;
buEdges(find(buEdges==0.0001))=0;
buEdges(find(buEdges>0.0001))=1;
degree(th)=mean(sum(buEdges));

thresh_P = [90 92.5 95 97.5];

%% active only

%density(th)=sum(sum(active_edges))/(length(active_edges)*(length(active_edges)-1));
%degree(th)=mean(sum(active_edges));

 end
 
subplot(1,2,2)
     scatter(thresholds, degree,'+','b');
    xlabel('STTC threshold') 
    ylabel('mean node degree')
    lineThickness(1.5) 
    aesthetics
    
