function Phase_Locking_Analysis(animals,KO,Frange,spthr,cols)

bnd = {'Delta';'Theta';'Beta';'Sgam';'Fgam'};
titOx = {'Pre-Oxy', 'Post-Oxy'};
titM = {'Motion', 'Immobility'};
titC = {'PY','IN'};

screensize = get(0,'Screensize');

la = length(animals);
lfr = size(Frange,1);

PhSpace = linspace(-90,270);

%% POOL UNIT TYPES FROM ALL SETS
load(fullfile('..','Analysis Results','Clusters.mat'));
for a = 1:la
    clustered_units{a} = cell2mat(clustered_units{a});
    for py = 1:2
        Nctype(a,py) = sum(clustered_units{a} == py);
    end
end

%% REMOVE BAD UNITS CONCATENATE UNITS FROM ALL SETS (ALONG 4TH DIMENSION) AND SPLIT CELL TYPES
mPhases = cell(la,2);
PhaseVec = cell(la,2);
PhaseRt = cell(la,2);

for a = 1:la
    S = load(fullfile('..','Analysis Results',animals{a},'Units.mat'));     % Load animal unit file
    S = remove_bad_units(animals{a},S);                                     % REMOVE PRE-SELECTED BAD UNITS
    
    S.mPhases = permute(S.mPhases,[2,3,4,1]);           % Stack sets (row cells) along the 4th dimension (each set: [f,x,m,units])
    S.PhaseVec = permute(S.PhaseVec,[2,3,4,1]);
    S.PhaseRt = permute(S.PhaseRt,[2,3,4,1]);
    
    S.mPhases = cell2mat(S.mPhases);                    % Concatenate sets (concatenates all units along 4th dimenion)
    S.PhaseVec = cell2mat(S.PhaseVec);
    S.PhaseRt = cell2mat(S.PhaseRt);
    
    for py = 1:2
        ct = find(clustered_units{a} == py);
        
        mPhases{a,py} = S.mPhases(:,:,:,ct);
        PhaseVec{a,py} = S.PhaseVec(:,:,:,ct);
        PhaseRt{a,py} = S.PhaseRt(:,:,:,ct);
    end
end
clear S

%% FDR-CORRECT RAYLEIGH PVALUES
% FDRPhaseRt = cell(la,2);
% for py = 1:2                                                                % For each cell group
%     Rt = PhaseRt(:,py);                                                     % Keep all Rayleigh pvalues of those cells
%     Rt = permute(Rt,[2,3,4,1]);            % Stack all animals along 4th dimension
%     Rt = cell2mat(Rt);                      % Concatenate animals (concatenates all units along 4th dimenion)
%     for x = 1:2                                                         % For each Oxy-condition
%         for m = 1:2                                                     % For each motion-condition
%             for f = 1:lfr                                               % For each frequency range
%                 [~,p] = fdr(squeeze(Rt(f,x,m,:)));
%                 
%                 for a = 1:la
%                     FDRPhaseRt{a,py}(f,x,m,:) = p(1:Nctype(a,py));
%                     p(1:Nctype(a,py)) = [];
%                 end
%             end
%         end
%     end
% end

%% BONFERRONI-CORRECT RAYLEIGH PVALUES
BONPhaseRt = cell(la,2);
for a = 1:la
    for py = 1:2
        BONPhaseRt{a,py} = PhaseRt{a,py}*sum(Nctype(:,py));
    end
end

%% FIND CELLS WITH NOT ENOUGH SPIKING
load(fullfile('..','Analysis Results','Rates_0.5.mat'),'totspikes');  % GET TOTAL NUMBER OF SPIKES PER UNIT PER CONDITION

lowspikers = cell(la,2,2,2);
for i = 1:numel(totspikes)
    lowspikers{i} = (totspikes{i} < spthr);
end

%% POOL PHASE LOCKING DATA (SPARSE SPIKING CONDITIONS AND NON-PHASE LOCKED UNITS)
mPh = cell(la,2,2,2,lfr);
PhV = cell(la,2,2,2,lfr);
for a = 1:la
    for py = 1:2                                                            % For each cell group
        for x = 1%:2                                                         % For each Oxy-condition
            for m = 1:2                                                     % For each motion-condition
                lsp = lowspikers{a,py,x,m};
                for f = 1:lfr                                               % For each frequency range
                    Rts = squeeze(BONPhaseRt{a,py}(f,x,m,:));
                    units = ~lsp' & Rts < 0.05;
                    
                    mPh{a,py,x,m,f} = [mPh{a,py,x,m,f}; squeeze(mPhases{a,py}(f,x,m,units))]; % Pool mean phases
                    PhV{a,py,x,m,f} = [PhV{a,py,x,m,f}; squeeze(PhaseVec{a,py}(f,x,m,units))];% and R-vectors

                    lockedcells(a,py,x,m,f) = sum(units);                       % Count number of cells to be included
                end
            end
        end
    end
end

%% POOL CELLS FROM ALL ANIMALS AND SPLIT GROUPS
MWT = cell(2,2,2,lfr);
MKO = cell(2,2,2,lfr);
VWT = cell(2,2,2,lfr);
VKO = cell(2,2,2,lfr);
for py = 1:2                                                                % For each cell group
    for x = 1%:2                                                             % For each oxy-condition
        for m = 1:2                                                         % For each motion-condition
            for f = 1:lfr                                                   % For each frequency range                
                MWT{py,x,m,f} = cell2mat(mPh(KO == 1,py,x,m,f));                                    % Keep mean phase locking of WTs
                MKO{py,x,m,f} = cell2mat(mPh(KO == 2,py,x,m,f));                                    % and KOs
                VWT{py,x,m,f} = cell2mat(PhV(KO == 1,py,x,m,f));                                    % Keep mean phase locking of WTs
                VKO{py,x,m,f} = cell2mat(PhV(KO == 2,py,x,m,f));                                    % and KOs
            end
        end
    end
end

%% COMPUTE MEANS PER ANIMAL AND SPLIT GROUPS
% MWT = cell(2,2,2,lfr);
% MKO = cell(2,2,2,lfr);
% VWT = cell(2,2,2,lfr);
% VKO = cell(2,2,2,lfr);
% for py = 1:2                                                                % For each cell group
%     for x = 1%:2                                                             % For each oxy-condition
%         for m = 1:2                                                         % For each motion-condition
%             for f = 1:lfr                                                   % For each frequency range                
%                 MWT{py,x,m,f} = cellfun(@mean,mPh(KO == 1,py,x,m,f));                                    % Keep mean phase locking of WTs
%                 MKO{py,x,m,f} = cellfun(@mean,mPh(KO == 2,py,x,m,f));                                    % and KOs
%                 VWT{py,x,m,f} = cellfun(@mean,PhV(KO == 1,py,x,m,f));                                    % Keep mean phase locking of WTs
%                 VKO{py,x,m,f} = cellfun(@mean,PhV(KO == 2,py,x,m,f));                                    % and KOs
%                 
%                 MWT{py,x,m,f}(isnan(MWT{py,x,m,f})) = [];    
%                 MKO{py,x,m,f}(isnan(MKO{py,x,m,f})) = [];                              
%                 VWT{py,x,m,f}(isnan(VWT{py,x,m,f})) = [];                                 
%                 VKO{py,x,m,f}(isnan(VKO{py,x,m,f})) = [];
%             end
%         end
%     end
% end

%% SIGNIFICANCE OF DISTRIBUTION PHASE LOCKING
pPhWTKO = ones(2,2,2,lfr);                                                  % Pvalues of mean phase locking of WTs vs KOs
pVWTKO = ones(2,2,2,lfr);                                                   % Pvalues of mean R-lengths of WTs vs KOs
pWT = ones(2,2,2,lfr);                                                      % Pvalues of individual phase locking of WT distributions over frequencies
pKO = ones(2,2,2,lfr);                                                      % Pvalues of individual phase locking of KO distributions over frequencies

for py = 1:2                                                                % For each cell group
    for x = 1%:2                                                             % For each oxy-condition
        for m = 1:2                                                         % For each motion-condition
            for f = 1:lfr                                                   % For each frequency range
                mWT = MWT{py,x,m,f};                                    % Keep mean phase locking of WTs
                mKO = MKO{py,x,m,f};                                    % and KOs
                vWT = VWT{py,x,m,f};                                    % Keep mean phase locking of WTs
                vKO = VKO{py,x,m,f};                                    % and KOs
                
                if ~isempty(mWT) && ~isempty(mKO)
                    pWT(py,x,m,f) = circ_rtest(mWT);                      % Rayleigh test for non-uniformity of the means of all WT units
                    pKO(py,x,m,f) = circ_rtest(mKO);                      % and KOs
                    
                    pWT(py,x,m,f) = pWT(py,x,m,f) * 2;                % BONFERRONI CORRECTION over the two distribution types
                    pKO(py,x,m,f) = pKO(py,x,m,f) * 2;                % (assuming each frequency has independent comparisons)
                    
                    if pWT(py,x,m,f) < 0.05 & pKO(py,x,m,f) < 0.05          % If both distributions are significantly phase locked on average
                        [pPhWTKO(py,x,m,f),~] = circ_wwtest(mWT,mKO); % Compare their mean phases with parametric Watson-Williams multi-sample test for equal means.
                    end
                end
                [pVWTKO(py,x,m,f),testtype] = significance(vWT,vKO,'unequal');testtype;  % Compare their mean R-values
            end
%             [~,p] = fdr(pPhWTKO(py,x,m,:));
%             pPhWTKO(py,x,m,:) = p;
%             
%             [~,p] = fdr(pVWTKO(py,x,m,:));
%             pVWTKO(py,x,m,:) = p;
        end
    end
end

% pPhWTKO=pPhWTKO*lfr;
% pVWTKO=pVWTKO*lfr;

%% COMPUTE DISTRIBUTION MEANS AND SEs
mP = zeros(2,2,2,2,lfr);
sP = zeros(2,2,2,2,lfr);
mV = zeros(2,2,2,2,lfr);
sV = zeros(2,2,2,2,lfr);
for py = 1:2                                                                % For each cell group
    for x = 1%:2                                                             % For each oxy-condition
        for m = 1:2                                                         % For each motion-condition
            for f = 1:lfr
                mWT = MWT{py,x,m,f};                                    % Keep mean phase locking of WTs
                mKO = MKO{py,x,m,f};                                    % and KOs
                vWT = VWT{py,x,m,f};
                vKO = VKO{py,x,m,f};
                
                mP(1,py,x,m,f) = circ_mean(mWT);                     % Get circular mean phase locking of all significantly locked cells
                [~,s] = circ_std(mWT);                      % and circular standard deviation
                sP(1,py,x,m,f) = s/sqrt(length(mWT));
                mP(2,py,x,m,f) = circ_mean(mKO);                     % Get circular mean phase locking of all significantly locked cells
                [~,s] = circ_std(mKO);                      % and circular standard deviation
                sP(2,py,x,m,f) = s/sqrt(length(mKO));
                
                mP(:,py,x,m,f) = rad2deg(mP(:,py,x,m,f));                                           % Turn to degrees
                sP(:,py,x,m,f) = rad2deg(sP(:,py,x,m,f));
                
                k = mP(:,py,x,m,f) < PhSpace(1);
                mP(k,py,x,m,f) = mP(k,py,x,m,f) + 360;
                
                mV(1,py,x,m,f) = mean(vWT);                          % and mean R-vector length
                sV(1,py,x,m,f) = std(vWT)/sqrt(length(vWT));
                mV(2,py,x,m,f) = mean(vKO);                          % and mean R-vector length
                sV(2,py,x,m,f) = std(vKO)/sqrt(length(vKO));
            end
        end
    end
end

%% PLOT PHASE LOCKING
count = 0;
for py = 1:2                                                                % For each cell group
    for x = 1%:2                                                             % For each oxy-condition
        for m = 1:2                                                         % For each motion-condition
            figure('Position', [30+count 30 screensize(3)/4 screensize(4)*0.7]);
            for f = 1:lfr
                subindex = (4*(4*(f-1)+1)+[1 2 3])+[0;4;8];
                subplot(4*lfr,4, subindex(:));hold on;
                title([num2str(sum(lockedcells(KO==1,py,x,m,f))), '  ',num2str(sum(lockedcells(KO==2,py,x,m,f)))]);
                
                mWT = MWT{py,x,m,f};                                    % Keep mean phase locking of WTs
                mKO = MKO{py,x,m,f};                                    % and KOs
                vWT = VWT{py,x,m,f};
                vKO = VKO{py,x,m,f};
                
                % Make cyclical by adding one more cycle
                PH1 = rad2deg(mWT);
                PH1 = [PH1; PH1+360];
                PH2 = rad2deg(mKO);
                PH2 = [PH2; PH2+360];
                
                % Plot distributions
                if py == 1
                    plot(PH1,repmat(vWT,2,1),'vb','Markersize',3);
                    plot(PH2,repmat(vKO,2,1),'vr','Markersize',3);
                else
                    plot(PH1,repmat(vWT,2,1),'ob','Markersize',3);
                    plot(PH2,repmat(vKO,2,1),'or','Markersize',3);
                end
                
                % Plot oscillation signal
                mx = max([vWT; vKO; 0]);                                    % Added 0 in case both PhV are empty
                plot(PhSpace, cos(deg2rad(PhSpace)+pi)*mx/2 + mx/2,'k','Linewidth',1.5);
                ylabel([bnd{f},'(',num2str(Frange(f,1)),'-',num2str(Frange(f,2)),'Hz)']);
                
                if f == lfr, xlabel('Phases');end;
                xlim([PhSpace(1) PhSpace(end)]);
                ylim([-0.00001 mx]);                                        % In case mx = 0
                box on;
                
                % Plot means and SEs of distributions
                if py == 1                                                  % For PYs
                    h1 = errorbar(mP(1,py,x,m,f),mV(1,py,x,m,f),...
                        sV(1,py,x,m,f),sV(1,py,x,m,f),sP(1,py,x,m,f),sP(1,py,x,m,f),...
                        's','Color','b','MarkerSize',10,'MarkerEdgeColor','b','Linewidth',2);
                    if pWT(py,x,m,f) < 0.05, set(h1,'MarkerFaceColor','b'); end
                    
                    h2 = errorbar(mP(2,py,x,m,f),mV(2,py,x,m,f),...
                        sV(2,py,x,m,f),sV(2,py,x,m,f),sP(2,py,x,m,f),sP(2,py,x,m,f),...
                        's','Color','r','MarkerSize',10,'MarkerEdgeColor','r','Linewidth',2);
                    if pKO(py,x,m,f) < 0.05, set(h2,'MarkerFaceColor','r'); end
                    
                else                                                        % For INs
                    h1 = errorbar(mP(1,py,x,m,f),mV(1,py,x,m,f),...
                        sV(1,py,x,m,f),sV(1,py,x,m,f),sP(1,py,x,m,f),sP(1,py,x,m,f),...
                        's','Color','b','MarkerSize',10,'MarkerEdgeColor','b','Linewidth',2);
                    if pWT(py,x,m,f) < 0.05, set(h1,'MarkerFaceColor','b'); end
                    
                    h2 = errorbar(mP(2,py,x,m,f),mV(2,py,x,m,f),...
                        sV(2,py,x,m,f),sV(2,py,x,m,f),sP(2,py,x,m,f),sP(2,py,x,m,f),...
                        's','Color','r','MarkerSize',10,'MarkerEdgeColor','r','Linewidth',2);
                    if pKO(py,x,m,f) < 0.05, set(h2,'MarkerFaceColor','r'); end
                end

                % Plot top histogram
                subplot(4*lfr,4, 4*4*(f-1)+[1 2 3]);hold on;
                h1 = histogram(PH1,PhSpace(1):10:PhSpace(end));
                h2 = histogram(PH2,PhSpace(1):10:PhSpace(end));
                Z = max([h1.Values, h2.Values]);
                line(mP(1,py,x,m,f)*[1 1],[0 Z], 'Color','b','Linewidth',2);
                line(mP(2,py,x,m,f)*[1 1],[0 Z], 'Color','r','Linewidth',2);
                plot_hist_signif(pPhWTKO(py,x,m,f),mP(1,py,x,m,f),mP(2,py,x,m,f),Z,Z,0,Z/10);
                axis tight;
                if f == 1, title([titC{py},'  ',titOx{x},'  ',titM{m},'  ','WT vs KO']); end
                set(gca,'Xtick',[]);
                
                % Plot side histogram
                subplot(4*lfr,4, 4*4*(f-1)+4*[2 3 4]);hold on; view(90, -90)
                if mx == 0, mx = 1; end                                     % So that the next line doesnt crash if mx = 0
                h1 = histogram(vWT,0:mx/20:mx);
                h2 = histogram(vKO,0:mx/20:mx);
                Z = max([h1.Values, h2.Values]);
                line(mV(1,py,x,m,f)*[1 1],[0 Z], 'Color','b','Linewidth',2);
                line(mV(2,py,x,m,f)*[1 1],[0 Z], 'Color','r','Linewidth',2);
                plot_hist_signif(pVWTKO(py,x,m,f),mV(1,py,x,m,f),mV(2,py,x,m,f),Z,Z,0,Z/10);
                set(gca,'Xtick',[]);
                axis tight;
            end
            count = count + 50;
            drawnow;
        end
    end
end

%% NUMBER OF PHASE LOCKED CELLS
for x = 1%:2
    figure;
    pv = ones(2,2);
    for m = 1:2
        for py = 1:2
            subplot(2,2,(m-1)*2+py);hold on;
            for f = 1:lfr
                C1 = lockedcells(KO == 1,py,x,m,f);                                      % Keep ratio of phase locked cells in each animal
                C2 = lockedcells(KO == 2,py,x,m,f);
                [pv(1,2), testtype] = significance(C1,C2,'unequal'); 
                
                plot_mean_SE([f-0.2,f+0.2],{C1,C2},pv,cols);
                
                disp(['-----',titOx{x},' ',titM{m},' ',titC{py},'-----']);
                disp(['Frequency = ',num2str(Frange(f,1)),':',num2str(Frange(f,2)),' Hz']);
                disp(['WT: ',num2str(mean(C1)),' +- ', num2str(std(C1))]);
                disp(['KO: ',num2str(mean(C2)),' +- ', num2str(std(C2))]);
                disp(['Pvalue = ',num2str(pv(1,2)),',  TestType = ',num2str(testtype)]);
            end
            set(gca,'Xtick',1:lfr,'Xticklabel',bnd);
            title([titOx{x},' ',titM{m},': ',titC{py}]);
            ylabel('Ratio of significantly locked cells per animal');
            count = count+2;
        end
    end
end
