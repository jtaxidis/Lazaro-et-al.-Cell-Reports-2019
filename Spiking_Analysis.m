function Spiking_Analysis(animals,KO)

%% POOL CELLS AND SPIKING PROPERTIES
load(fullfile('..','Analysis Results','Clusters.mat'));

la = length(animals);

mR = cell(la,2,2,2);
BI = cell(la,2,2,2);
CSI = cell(la,2,2,2);
ISI = cell(la,2,2,2);

for a = 1:la                                                                % For each animal
    S = load(fullfile('..','Analysis Results',animals{a},'Units.mat'));     % Load animal unit file
    
    S = remove_bad_units(animals{a},S);                                     % REMOVE PRE-SELECTED BAD UNITS
    
    ls = length(S.mRate);                                                   % Number of sets for that animal
    
    for x = 1:2                                                             % For each condition
        for py = 1:2                                                        % For each cell type
            for m = 1:2                                                     % For motion/immobility
                mR{a,py,x,m} = [];
                BI{a,py,x,m} = [];
                CSI{a,py,x,m} = [];
                ISI{a,py,x,m} = [];
                
                for st = 1:ls                                               % For each set
                    cells = (clustered_units{a}{st} == py);                 % Get indexes of that type of cells
                    lc = sum(cells);
                    
                    k = squeeze(S.mRate{st}(x,m,cells));                    % Keep the firing rates in each codnition for that cell group
                    mR{a,py,x,m} = cat(1,mR{a,py,x,m}, k);                  % Concatenate them
                    
                    k = squeeze(S.BIndex{st}(x,m,cells));                   % Keep the burst indexes for the cells that group
                    BI{a,py,x,m} = cat(1,BI{a,py,x,m}, k);                  % Concatenate them over all sets
                    
                    k = squeeze(S.CSIndex{st}(x,m,cells));                  % Keep the CSI indexes for the cells that group
                    CSI{a,py,x,m}= cat(1,CSI{a,py,x,m}, k);                 % Concatenate them over all sets
                    
                    k = squeeze(S.InterSpike{st}(x,m,cells));               % Keep the interspike intervals for the cells that group
                    ISI{a,py,x,m}= cat(1,ISI{a,py,x,m}, k);                 % Concatenate them over all sets
                end
            end
        end
    end
end

%% PLOT COMPARISONS
cols = ['b','r'];
% plot_hist_comp(mR,KO,'Rates',cols)                                          % BAD VERSION OF MEAN RATE (REPLACED WITH MORE ACCURATE IN RATES_ANALYSIS)
% plot_hist_comp(ISI,KO,'ISI',cols)                                           % BAD VERSION OF ISI (REPLACED WITH MORE ACCURATE IN RATES_ANALYSIS)
plot_hist_comp(BI,KO,'BI',cols)
% plot_hist_comp(CSI,KO,'CSI',cols)

 