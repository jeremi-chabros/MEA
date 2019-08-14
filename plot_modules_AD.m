%plotting the clusters 

%load
clear adjM
load('MPT190515_2B_DIV28_cSpikes_L0_adjM.mat')
load('MPT190515_2B_DIV28_cSpikes_L0.mat')


%% whole array; weighted (below i could do this only on active adjM; and/or binary adjM)

%remove self connections
%adjM=adjM-eye(size(adjM));
%remove NaNs
adjM(isnan(adjM))=0;
%remove -ve weights
adjM(find(adjM<0))=0;
%set self connections to one
adjM=adjM+eye(size(adjM));

%plot
figure
imagesc(adjM)
%check modularity and get communities
[C,Q]=modularity_und(adjM);

%prepare for plot
[X,Y,INDSORT] = grid_communities(C); 
%remove self edge  2s
adjM(find(adjM>1))=1;

figure
imagesc(adjM(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);             % plot community boundaries
xticks([])
xticklabels([])
yticks([])
yticklabels([])
ylabel('Channels')
cb=colorbar
set(gca,'fontsize',24)

%plot only good clusters 1 2 38 39
figure
%only chosen clusters - ones that have more than just one self connection
imagesc(adjM(INDSORT(find(sum(adjM(INDSORT,INDSORT))~=1)),INDSORT(find(sum(adjM(INDSORT,INDSORT))~=1))))

%wanted communities=[1 2 38 39];
Cw=C(find(C==1 | C==2 | C==38 | C==39));
Aw=adjM(find(C==1 | C==2 | C==38 | C==39),find(C==1 | C==2 | C==38 | C==39));
figure 
imagesc(Aw)
[X,Y,INDSORT] = grid_communities(Cw);


figure
imagesc(Aw(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);             % plot community boundaries
xticks([])
xticklabels([])
yticks([])
yticklabels([])
ylabel('Channels')
cb=colorbar
set(gca,'fontsize',24)

%could compare Q to that in a matched random graph - need to find a WU
%graph randomising algorhithm 
[R,eff]=randmio_und(Aw, 100); %randomised weighted network
[Cr,Qr]=modularity_und(R); %randomisede parameters
[X,Y,INDSORT] = grid_communities(Cr); 

figure
imagesc(R(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);             % plot community boundaries
xticks([])
xticklabels([])
yticks([])
yticklabels([])
ylabel('Channels')
cb=colorbar
set(gca,'fontsize',24)


% GRID_COMMUNITIES       Outline communities along diagonal
%
%   [X Y INDSORT] = GRID_COMMUNITIES(C) takes a vector of community
%   assignments C and returns three output arguments for visualizing the
%   communities. The third is INDSORT, which is an ordering of the vertices
%   so that nodes with the same community assignment are next to one
%   another. The first two arguments are vectors that, when overlaid on the
%   adjacency matrix using the PLOT function, highlight the communities.
%
%   Example:
%
%   >> load AIJ;                                % load adjacency matrix
%   >> [C,Q] = modularity_louvain_und(AIJ);     % get community assignments
%   >> [X,Y,INDSORT] = fcn_grid_communities(C); % call function
%   >> imagesc(AIJ(INDSORT,INDSORT));           % plot ordered adjacency matrix
%   >> hold on;                                 % hold on to overlay community visualization
%   >> plot(X,Y,'r','linewidth',2);             % plot community boundaries
%
%   Inputs:     C,       community assignments
%
%   Outputs:    X,       x coor
%               Y,       y coor
%               INDSORT, indices
%
%   Richard Betzel, Indiana University, 2012

%% active only

clear adjM
load('MPT190515_2B_DIV28_cSpikes_L0_adjM.mat')
load('MPT190515_2B_DIV28_cSpikes_L0.mat')
%get active adjM
ac_chan_vec=full(find(sum(cSpikes)>9));
adjM=adjM(ac_chan_vec,ac_chan_vec);

%% whole array; weighted (below i could do this only on active adjM; and/or binary adjM)

%remove self connections
    adjM=adjM-eye(size(adjM));
%remove NaNs
adjM(isnan(adjM))=0;
%remove -ve weights
adjM(find(adjM<0))=0;
%set self connections to one
    %adjM=adjM+eye(size(adjM));

%plot
figure
imagesc(adjM)
%check modularity and get communities
gamma=1
[C,Q]=modularity_und(adjM,gamma);

%prepare for plot
[X,Y,INDSORT] = grid_communities(C); 

figure
imagesc(adjM(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);             % plot community boundaries
xticks([])
xticklabels([])
yticks([])
yticklabels([])
ylabel('Channels')
cb=colorbar
set(gca,'fontsize',24)