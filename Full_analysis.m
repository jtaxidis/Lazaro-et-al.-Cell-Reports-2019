%% PARAMETERS
addpath('CircStat2012a');

animals={'M2420';
    'M2422';
    'M2453';
    'M2455';
    'M2730';
    'M2870';                                                                % Weirdly low post-oxy LFP power
    
    'MD1A';
    'MD1B';
    'MD4B';
    'MD4C';
    'MD2A';
    'MD2B';
    'MD3A'};

KO = [2,1,1,1,2,1, 1,1,1,1,2,2,2];                                          % 1 = WT, 2 = KO

channels = [51 52 52 14 52 52   48,49,112 97 53 48 50];                     % Contains layers 2/3, 5, 6
% channels = [118 114 118 75 52 52   16,49,112 97 118 118 118];

bnd = {'Delta';'Theta';'Beta';'Sgam';'Fgam'};
Frange = [1,4;
    5,11;
    12,30;
    30,55;
    80,110];

titOx = {'Pre-Oxy', 'Post-Oxy'};
titM = {'Motion', 'Immobility'};
titC = {'PY','IN'};
titKO = {'WT','KO'};

spthr = 200;                                                                % MINIMUM NUMBER OF SPIKES PER CONDITION

cols = ['b';'r'];

%% GET MOTION SEGMENTS
% for m = 1:length(animals)
%     disp(['PERFORMING MOTION ANALYSIS ON ANIMAL: ',animals{m}]);
%     Get_Motion(animals{m});
%     disp(' ');
% end

%% UNIT ANALYSIS
% for m = 1:length(animals)           
%     disp(['PERFORMING UNIT ANALYSIS ON ANIMAL: ',animals{m}]);
%     Unit_Analysis(animals{m});
%     disp(' ');
% end

%% UNIT CLUSTERING (BAD UNITS ARE REMOVED HERE)
% Unit_Clustering(animals);
% drawnow;

%% COMPUTE FIRING RATES, ISI AND CORRELATIONS
binlen(1) = 0.5;    %delta 2Hz                                              % Time length of each rate bin
binlen(2) = 0.1;    %theta 10Hz                                                           
binlen(3) = 0.05;    %beta 20Hz                                                           
binlen(4) = 0.025;    %slow gamma 40Hz                                       
binlen(5) = 0.01; % gamma 100 Hz

for i = 1:5
    Get_Rates(animals,binlen(i),spthr);
end

%% MOTION ANALYSIS
Motion_Analysis(animals,KO);

%% SPIKING PROPERTIES ANALYSIS
Spiking_Analysis(animals,KO);

%% FIRING RATES ANALYSIS
Rates_Analysis(animals,KO,binlen(2));       % USE THIS FOR RATES/ISI/FANO

%% CORRELATIONS ANALYSIS
Rates_Analysis(animals,KO,binlen(1),[-0.4 0.9; -0.02, 1]);       % USE THIS FOR CORRELATIONS
Rates_Analysis(animals,KO,binlen(2),[-0.15 0.6; -0.02, 0.8]);      
Rates_Analysis(animals,KO,binlen(3),[-0.1 0.4; -0.02, 0.8]);
Rates_Analysis(animals,KO,binlen(4),[-0.05 0.4; -0.02, 0.55]);
Rates_Analysis(animals,KO,binlen(5),[-0.05 0.1; -0.02, 0.3]);

%% PHASE LOCKING ANALYSIS
Phase_Locking_Analysis(animals,KO,Frange,spthr,cols)


%% ------------------------- LFP ANALYSIS ---------------------------------
%% COMPUTE LFP BANDPOWER
for m = 1%:length(animals)
    disp(['COMPUTING LFP BANDPOWER ON ANIMAL: ',animals{m}]);
    LFP_bandpower(animals{m},channels(m),Frange);
    disp(' ');
end

%% POOL LFP BANDPOWER AND COMPUTE SIGNIFICANCE
lfr = size(Frange,1);
LFPsp = cell(2,2,2);
LFPpower = cell(2,2,2);
pval = cell(2,2);

for a = 1:length(animals)                                                   % For each animal
    S = load(fullfile('..','Analysis Results',animals{a},['LFPpower_ch',num2str(channels(a)),'.mat']));     % Load animal LFPpower file for that channel
    for x = 1:2
        for m = 1:2
            ind = S.indexes{x,m};                                           % Load the corresponding LFP indexes
            if sum(ind)*(1/1000) >= 30                                      % IF THE ANIMAL WAS AT THAT MOTION/IMMOBITILY STATE FOR AT LEAST A MINUTE IN TOTAL (to allow fo sufficient data)
                P = squeeze(S.LFPpower(x,m,:))';                                % Keep the corresponding bandpower over all frequencies
                P = log(P);                                                     % Keep log power (OPTIONAL)
                %             P = P / P(1);                                                   % Divide by Delta power (OPTIONAL)
                LFPpower{KO(a),x,m} = [LFPpower{KO(a),x,m}; P];                 % Stack as row (yields KO x lfr matrix)
                
                LFPsp{KO(a),x,m} = [LFPsp{KO(a),x,m}; S.LFPspec{x,m}'];         % Stack power spectra
            end
        end
    end
end
freq = S.freq;

for x = 1:2
    for m = 1:2
        pval{x,m} = ones(lfr,1);
        for f = 1:lfr
            [pval{x,m}(f),testtype] = significance(LFPpower{1,x,m}(:,f),LFPpower{2,x,m}(:,f),'unequal');testtype;
             pval{x,m}(f) = ranksum(LFPpower{1,x,m}(:,f),LFPpower{2,x,m}(:,f),'alpha',0.05,'tail','both');
        end
        pval{x,m} = pval{x,m} * lfr;                                        % BONFERRONI CORRECTION over all frequency comparions
    end
end

%% PLOT LFP BANDPOWER
figure;
count = 2;
for x = 1:2
    for m = 1:2
        subplot(4,2,count);hold on;
        for f = 1:lfr
%             plot_mean_SE([f-0.2,f+0.2],{LFPpower{1,x,m}(:,f),LFPpower{2,x,m}(:,f)},[1,pval{x,m}(f,:)],cols);
            plot(f*ones(size(LFPpower{1,x,m},1),1) , LFPpower{1,x,m}(:,f),'o','Color',cols(1));
            plot(f*ones(size(LFPpower{2,x,m},1),1) , LFPpower{2,x,m}(:,f),'o','Color',cols(2));
        end
        errorbar(1:lfr,mean(LFPpower{1,x,m},1),std(LFPpower{1,x,m},1)/sqrt(size(LFPpower{1,x,m},1)),'d-','Color',cols(1),'Linewidth',2);
        errorbar(1:lfr,mean(LFPpower{2,x,m},1),std(LFPpower{2,x,m},1)/sqrt(size(LFPpower{2,x,m},1)),'d-','Color',cols(2),'Linewidth',2);
        plot(find(pval{x,m} < 0.05), ones(sum(pval{x,m}<0.05),1),'*k');
        
        set(gca,'Xtick',1:lfr,'Xticklabel',bnd);
        title([titOx{x},' ',titM{m}]);
        ylabel('LFP bandpower (Log)');
        xlim([0.8 lfr+0.2]);
        count = count+2;
    end
end

%% PLOT LFP SPECTRA
count = 1;
for x = 1:2
    for m = 1:2
        subplot(4,2,count);hold on;
        for k = 1:2
            mLFP = mean(LFPsp{k,x,m},1);                                    % Get mean power spectrum
            p = polyfit(log(freq'),mLFP,1);                                 % Fit it linearly in a log-log scale
            LFPsp{k,x,m} = LFPsp{k,x,m} - repmat(p(1)*log(freq')+p(2), size(LFPsp{k,x,m},1),1); % Remove the linear trend from original LFP
            
            mLFP = mean(LFPsp{k,x,m},1);                                    % Recompute mean
            sLFP = std(LFPsp{k,x,m},[],1)/sqrt(size(LFPsp{k,x,m},1));       % And SE
            
            fill_plot(freq',mLFP,0,sLFP,cols(k));
            plot(freq, mLFP,'Linewidth',1.5,'Color',cols(k));
        end
        title([titOx{x},' ',titM{m}]);
        ylabel('LFP power (log)');
        xlabel('Frequency (Hz)');
        xlim([0 120]);
        count = count+2;
    end
end
