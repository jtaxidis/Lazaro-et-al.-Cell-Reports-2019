function Unit_Clustering(animals)

Fs = 25e3;
window = 35;
time = (-window:window)/Fs * 1000;   %time axis in ms

titM = {'Motion', 'Immobility'};

Na = length(animals);

ls = zeros(Na,1);                                       
Nu = cell(Na,1);                

rng(1);                                                                     % For reproducibility

%% ALLOCATE MEMORY
R = cell(2,2);
BI = cell(2,2);
CSI = cell(2,2);
for x = 1:2
    for m = 1:2
        R{x,m} = [];
        BI{x,m} = [];
        CSI{x,m} = [];
    end
end
PT = [];
PTr = [];
WF = [];
HW = [];

%% POOL ANIMALS
for a = 1:Na
    S = load(fullfile('..','Analysis Results',animals{a},'Units.mat'));     % Load animal unit file   
    S = remove_bad_units(animals{a},S);                                     % REMOVE PRE-SELECTED BAD UNITS

    Peak = S.Peak;                                                          % Keap all relevant measures
    Trough = S.Trough;
    PTdist = S.PTdist;
    Waveforms = S.Waveforms;
    Halfwidth = S.Halfwidth;
    mRate = S.mRate;
    BIndex = S.BIndex;
    CSIndex = S.CSIndex;
    
    ls(a) = length(Peak);                                                   % Number of sets
    Nu{a} = cellfun(@length,Peak);                                          % number of units in each set
    
    for x = 1:2                                                             % For each oxy-condition
        for m = 1:2                                                         % For each motion-condition
            for st = 1:ls(a)                                                % For each set
                R{x,m} = [R{x,m}; squeeze(mRate{st}(x,m,:))];               % Keep corresponding mean rates
                BI{x,m} = [BI{x,m}; squeeze(BIndex{st}(x,m,:))];            % BIs
                BI{x,m}(isnan(BI{x,m})) = 0;
                CSI{x,m} = [CSI{x,m}; squeeze(CSIndex{st}(x,m,:))];         % CSIs
                CSI{x,m}(isnan(CSI{x,m})) = 0;
            end
        end
    end
    PT = [PT;cell2mat(PTdist)];                                             % And all Peak-Trough distances
    PTr = [PTr; cell2mat(Peak)./abs(cell2mat(Trough))];                     % Peak-trough ratios
    WF = [WF; cell2mat(Waveforms)]; 
    HW = [HW;cell2mat(Halfwidth)];
end

%% DIFFERENT CLUSTERING METHODS (UNSUCCESSFUL CLUSTERING)
% X = [PT,PTr];
% figure;

% %% K-MEANS
% [idx,C] = kmeans(X(:,1:2),2,'Replicates',10,'Start','plus');
% [~,smallclust] = min([sum(idx==1),sum(idx==2)]);
% bigclust = setdiff(1:2,smallclust);
% C1 = find(idx == bigclust);
% C2 = find(idx == smallclust);
% plot_clusters(X,WF,C1,C2,time,[],1,'K-means');
% 
% %% SOM
% addpath('SOM/somtoolbox');
% D = som_data_struct(WF);                       % Turn LFP patterns to appropriate data structure                                                                       
% D = som_normalize(D,'var');                                                 % Normalize data
% sM = som_make(D,'msize',2,'lattice','rect','training','long');
% sM = som_autolabel(sM,D,'freq');                                            % Add labels to SOM
% bmu = som_bmus(sM,D);                                                       % Best matching unit for each LFP pattern
% [~,smallclust] = min([sum(bmu==1),sum(bmu==2)]);
% bigclust = setdiff(1:2,smallclust);
% C1 = find(bmu == bigclust);
% C2 = find(bmu == smallclust);
% plot_clusters(X,WF,C1,C2,time,[],2,'Self Organizing Map');
% 
% % hits = som_hits(sM,D);                                                      % Number of hits in each map unit
% % units = sM.codebook;                                                        % Coordinates of each unit
% % bmu_coord = units(bmu,:);                                                   % Coordinates of each pattern's BMU
% 
% %% DENDROGRAM
% WFz = zscore(WF);
% Z = linkage(WFz,'complete','euclidean');         % Because of euclidean distance 1-zLFP_corrs is the same (due to squaring differences)
% idx = cluster(Z,'maxclust',2,'criterion','distance');
% % dendrogram(Z,0,'colorthreshold',30);
% [~,smallclust] = min([sum(idx==1),sum(idx==2)]);
% bigclust = setdiff(1:2,smallclust);
% C1 = find(idx == bigclust);
% C2 = find(idx == smallclust);
% plot_clusters(X,WF,C1,C2,time,[],3,'Hierarchical cluster tree');

% %% PCA & KMEANS
% Y = [PT,PTr,BI{1,1},CSI{1,1}];              % Combine Peak-trough distance, ratio, burst index and SCI 
% Yz = zeros(size(Y));
% for c = 1:size(Y,2)
%     Yz(:,c) = zscore(Y(:,c));                % z-score
% end
% [~,score,~,~,explained] = pca(Yz);           % PCA 
% Ypca = score(:,1:3);                        % Keep first 2 PCs
% expl = explained(1:2);
% [idx,C] = kmeans(Ypca(:,1:2),2,'Replicates',100,'Start','plus');  % Cluster them with kmeans
% [~,smallclust] = min([sum(idx==1),sum(idx==2)]);
% bigclust = setdiff(1:2,smallclust);
% C1 = find(idx == bigclust);
% C2 = find(idx == smallclust);
% idx(C1) = 1;                % SET THE LARGER CLUSTER AS CLUSTER 1
% idx(C2) = 2;
% plot_clusters(X,WF,C1,C2,time,[],4,'PCA & K-means');

%% PCA & KMEANS
Y = [PT,PTr,BI{1,1}];                                                       % Combine Peak-trough distance, Peak-trough ratio and pre-oxy motion burst index 
Yz = zscore(Y,[],1);                                                        % z-score each variable

[~,score,~,~,explained] = pca(Yz);                                          % PCA 
Ypca = score(:,1:3);                                                        % Keep first 3 PCs (ALL OF THEM)
expl = explained(1:3)
[idx,C] = kmeans(Ypca,2,'Replicates',100,'Start','plus');                   % Cluster PCs into 2 cluster with kmeans

[~,smallclust] = min([sum(idx==1),sum(idx==2)]);                            % Find the smallest cluster
bigclust = setdiff(1:2,smallclust);
C1 = find(idx == bigclust);
C2 = find(idx == smallclust);
idx(C1) = 1;                                                                % SET THE LARGER CLUSTER AS CLUSTER 1
idx(C2) = 2;

%% MAP BACK TO INDIVIDUAL ANIMALS
idxtemp = idx;                                                              % Reproduce the cluster indexes (pooled all animals, sets, units)
clustered_units = cell(Na,1);       
for a = 1:Na                                                                % For each animal
    clustered_units{a} = cell(ls(a),1); 
    for st = 1:ls(a)                                                        % For each set
        clustered_units{a}{st} = idxtemp(1:Nu{a}(st));                      % Store the cluster indexes of the first Nu units (equal to the units of that set)
                                                                            % (same order as when doing cell2mat and pooling)
        idxtemp(1:Nu{a}(st)) = [];                                          % Delete the corresponding indexes
    end
end

%% PLOT PCA IN MORE DETAIL
figure;
plot_clusters(Y(:,1:2),WF,C1,C2,time,[],1,'PCA & K-means');

subplot(2,4,2); hold on;
plot3(Y(C1,1),Y(C1,2),Y(C1,3),'ok','markerfacecolor','k');
plot3(Y(C2,1),Y(C2,2),Y(C2,3),'or','markerfacecolor','r'); view(3); axis square; grid on;
title('PCA & K-means'); xlabel('PT-dist');ylabel('PT-ratio');zlabel('BI');

% subplot(3,4,[3 4 7 8]); hold on;
% plot3(Y(C1,1),Y(C1,2),Y(C1,4),'ok','markerfacecolor','k');
% plot3(Y(C2,1),Y(C2,2),Y(C2,4),'or','markerfacecolor','r'); view(3); axis square; grid on;
% title('PCA & K-means'); xlabel('PT-dist');ylabel('PT-ratio');zlabel('CSI');

subplot(2,4,3);hold on;
plot3(Ypca(C1,1),Ypca(C1,2),Ypca(C1,3),'ok','markerfacecolor','k');
plot3(Ypca(C2,1),Ypca(C2,2),Ypca(C2,3),'or','markerfacecolor','r');  view(3); axis square; grid on;
title('PCA space'); xlabel('PCA1');ylabel('PCA2');  zlabel('PCA3');

subplot(2,4,6); hold on;
compare_hist(PT,C1,C2);
xlabel('Peak-Trough distance (ms)'); 

subplot(2,4,7); hold on;
compare_hist(PTr,C1,C2);
xlabel('Peak-Trough Ratio'); axis tight;

subplot(2,4,8); hold on;
compare_hist(HW,C1,C2);
xlabel('Half Width (ms)'); axis tight;

%% PLOT FIRING RATE PROPERTIES
figure;
subplot(331); hold on;
plot(R{1,1}(C1),R{1,2}(C1),'ok','Markerfacecolor','k');
plot(R{1,1}(C2),R{1,2}(C2),'or','Markerfacecolor','r');
maxR = max([R{1,1};R{1,2}]);
line([0 maxR],[0 maxR],'Color','k');
title('Firing Rates');
xlabel(titM{1});ylabel(titM{2}); axis square; axis tight;

subplot(332); hold on;
plot(BI{1,1}(C1),BI{1,2}(C1),'ok','Markerfacecolor','k');
plot(BI{1,1}(C2),BI{1,2}(C2),'or','Markerfacecolor','r');
maxBI = max([BI{1,1};BI{1,2}]);
line([0 maxBI],[0 maxBI],'Color','k');
title('Burst Index');
xlabel(titM{1});ylabel(titM{2}); axis square; axis tight;

subplot(333); hold on;
plot(CSI{1,1}(C1),CSI{1,2}(C1),'ok','Markerfacecolor','k');   
plot(CSI{1,1}(C2),CSI{1,2}(C2),'or','Markerfacecolor','r');  
maxC = max([CSI{1,1};CSI{1,2}]);
minC = min([CSI{1,1};CSI{1,2}]);
line([minC maxC],[minC maxC],'Color','k');
title('SCI');
xlabel(titM{1});ylabel(titM{2}); axis square; axis tight;

disp(['Mean BI in two clusters during motion = ', num2str(mean(BI{1,1}(C1))),', ', num2str(mean(BI{1,1}(C2)))]); 
disp(['Mean BI in two clusters during immobility = ', num2str(mean(BI{1,2}(C1))),', ', num2str(mean(BI{1,2}(C2)))]); 

%% PLOT HISTOGRAMS
subplot(334); hold on;
compare_hist(R{1,1},C1,C2);
title('Firing Rates during Motion (pre-Oxy)'); 

subplot(337); hold on;
compare_hist(R{1,2},C1,C2);
title('Firing Rates during Immobility (pre-Oxy)');

subplot(335); hold on;
compare_hist(BI{1,1},C1,C2);
title('Burst Index during Motion (pre-Oxy)');

subplot(338); hold on;
compare_hist(BI{1,2},C1,C2);
title('Burst Index during Immobility (pre-Oxy)');

subplot(336); hold on;
compare_hist(CSI{1,1},C1,C2);
title('CSI during Motion (pre-Oxy)');

subplot(339); hold on;
compare_hist(CSI{1,2},C1,C2);
title('CSI during Immobility (pre-Oxy)');

%% SAVE
% save(fullfile('..','Analysis Results','Clusters.mat'),'clustered_units',...
%        'Y','Ypca','expl','idx','C','C1','C2');
   
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function plot_clusters(X,WF,C1,C2,time,C,sub,tit)
subplot(2,4,sub); hold on;
plot(X(C1,1),X(C1,2),'o','Color','k','Markerfacecolor','k');                % Plot the large cluster first
plot(X(C2,1),X(C2,2),'o','Color','r','Markerfacecolor','r');
if ~isempty(C)
    plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',2);
end
xlabel('Peak-Trough distance');
ylabel('Peak-Trough ratio');
title(tit);

subplot(2,4,sub+4); hold on;
fill_plot(time,mean(WF(C1,:)),0,std(WF(C1,:)),'k')
plot(time,mean(WF(C1,:)),'k','linewidth',2);
fill_plot(time,mean(WF(C2,:)),0,std(WF(C2,:)),'r')
plot(time,mean(WF(C2,:)),'r','linewidth',2);
xlabel('time (ms)');
ylabel('Average waveforms');
lc1 = length(C1);
lc2 = length(C2);
title(['PYs = ',num2str(lc1),' (',num2str(100*lc1/(lc1+lc2)),'%). INs = ',num2str(lc2),' (',num2str(100*lc2/(lc1+lc2)),'%']);
axis tight

function compare_hist(X,C1,C2)
X1 = X(C1);
X2 = X(C2);
m1 = min(X)-0.1;
m2 = max(X)+0.1;
binstep = (m2-m1)/25;                                                       % 30 bins
histogram(X1,m1:binstep:m2,'Facecolor','k');   
histogram(X2,m1:binstep:m2,'Facecolor','r');   
axis tight;


