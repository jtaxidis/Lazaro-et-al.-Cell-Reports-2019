function Rates_Analysis(animals,KO,binlen,corrlims)

load(fullfile('..','Analysis Results',['Rates_',num2str(binlen),'.mat']));

la = length(animals);

titOx = {'Pre-Oxy', 'Post-Oxy'};
titM = {'Motion', 'Immobility'};
titC = {'PY','IN'};
titKO = {'WT','KO'};

clim = [-0.5 0.5];

cols = ['b','r'];

%% PLOT FIRING RATE COMPARISON
plot_hist_comp(mRate,KO,'Rates',cols)

%% PLOT ISI COMPARISON
% plot_hist_comp(ISI,KO,'ISI',cols)

%% PLOT FANO FACTOR COMPARISON
FF = {};
for a = 1:la
    for x = 1%:2
        for m = 1:2
            for py = 1:2
                FF{a,py,x,m} = sRate{a,py,x,m}./mRate{a,py,x,m};
            end
        end
    end
end
% plot_hist_comp(FF,KO,'Fano Factor',cols)


% %% PLOT CORRELATION MATRICES
% for x = 1%:2
%     for m = 1:2
%         figure;
%         for py = 1:2
%             Rwt = Rho(KO==1,py,x,m);                                       % Keep WT mice
%             Rko = Rho(KO==2,py,x,m);                                       % And KO mice
%             
%             for a = 1:length(Rwt)
%                 Rwt{a}(isnan(Rwt{a})) = 0;
%             end
%             for a = 1:length(Rko)
%                 Rko{a}(isnan(Rko{a})) = 0;
%             end
%             
%             if m == 1
%                 clim = [-0.5 0.9];
%             else
%                 clim = [-0.02 1];
%             end
%             subplot(3,2,(py-1)*2+1)
%             imagesc(blkdiag(Rwt{1},Rwt{2},Rwt{3},Rwt{4},Rwt{5},Rwt{6},Rwt{7},Rwt{8}),clim);
%             axis square;
%             title(['WT ',titM{m},': ',titC{py}]);
%             
%             subplot(3,2,(py-1)*2+2)
%             imagesc(blkdiag(Rko{1},Rko{2},Rko{3},Rko{4},Rko{5}),clim);
%             axis square;
%             title(['KO ',titM{m},': ',titC{py}]);
%         end
%         
%         cRwt = cRho(KO==1,x,m);                                       % Keep WT mice
%         cRko = cRho(KO==2,x,m);                                       % Keep WT mice
%         
%         for a = 1:length(cRwt)
%             cRwt{a}(isnan(cRwt{a})) = 0;
%         end
%         for a = 1:length(cRko)
%             cRko{a}(isnan(cRko{a})) = 0;
%         end
%         
%         subplot(3,2,5)
%         imagesc(blkdiag(cRwt{1},cRwt{2},cRwt{3},cRwt{4},cRwt{5},cRwt{6},cRwt{7},cRwt{8})',clim);
%         axis square;
%         title(['WT ',titM{m},': PY-IN']);
%         
%         subplot(3,2,6)
%         imagesc(blkdiag(cRko{1},cRko{2},cRko{3},cRko{4},cRko{5})',clim);
%         axis square;
%         title(['KO ',titM{m},':PY-IN']);
%     end
% end
%

%% PLOT NUMBER OF CELLS THAT SPIKED ENOUGH TO BE INCLUDED IN CORRELATIONS
cells = zeros(la,2,2,2);
for a = 1:length(animals)
    for x = 1%:2
        for m = 1:2
            for py = 1:2
                A = diag(Rho{a,py,x,m},1);
                A(isnan(A)) = [];
                cells(a,py,x,m) = length(A);%/size(Rho{a,py,x,m},1) * 100;                 
            end
        end
    end
end

for x = 1%:2
    figure('Name','Correlations');   
    for m = 1:2                                                             % For motion/immobility
        for py = 1:2                                                        % For each cell type
            subplot(2,2,(m-1)*2+py)
            
            [pvalue, testtype] = significance(cells(KO == 1,py,x,m), cells(KO == 2,py,x,m),'unequal'); 
            pv = ones(2,2);
            pv(1,2) = pvalue;% * 2;                                           % BONFERRONI CORRECTION OVER 2 CELL TYPE COMPARISONS

            plot_mean_SE(1:2,{cells(KO == 1,py,x,m), cells(KO == 2,py,x,m)},pv,cols);
                       
            disp(['-----',titOx{x},' ',titM{m},' ',titC{py},'-----']);
            disp([titKO{1},': ',num2str(mean(cells(KO == 1,py,x,m))),' +- ', num2str(std(cells(KO == 1,py,x,m)))]);
            disp([titKO{2},': ',num2str(mean(cells(KO == 2,py,x,m))),' +- ', num2str(std(cells(KO == 2,py,x,m)))]);
            disp(['Pvalue = ',num2str(pv(1,2)),',  TestType = ',num2str(testtype)]);
        end
    end
end

            
            
%% PLOT CORRELATION DISTRIBUTIONS
for a = 1:length(animals)
    for x = 1%:2
        for m = 1:2
            for py = 1:2
                Rho{a,py,x,m} = full(upper_triangular(Rho{a,py,x,m}));
                Rho{a,py,x,m}(isnan(Rho{a,py,x,m})) = [];                   % Remove correlations between cells that never spiked or didnt spike enough
            end
            cRho{a,x,m} = full(upper_triangular(cRho{a,x,m}));
            cRho{a,x,m}(isnan(cRho{a,x,m})) = [];
        end
    end
end

nbins = 50;                                                                 % Number of histogram bins
for x = 1%:2
    figure('Name','Correlations');   
    for m = 1:2                                                             % For motion/immobility
        for py = 1:2                                                        % For each cell type
            rwt = Rho(KO==1,py,x,m);                                          % Keep WT mice
            rko = Rho(KO==2,py,x,m);                                          % And KO mice
            
            if size(rwt{8},2) == 0 % Fix one empty problematic correlation matrix
                rwt{8} = rwt{8}';
            end
            
            % POOL ALL CELLS ------------------------------------------
            RWT = cell2mat(rwt);                                            % Concatenate cells from all animals
            RKO = cell2mat(rko);
            
            RWT(isnan(RWT)) = [];
            RKO(isnan(RKO)) = [];
            
            [pvalue, testtype] = significance(RWT,RKO,'unequal'); testtype;
            pvalue = pvalue * 3;                                            % BONFERRONI CORRECTION OVER 3 CELL TYPE COMPARISONS
            
            subplot(2,9,(m-1)*9+(py-1)*3+[1 2]);hold on;
            m1 = corrlims(m,1);
            m2 = corrlims(m,2);
            binstep = (m2-m1)/nbins;                                       
            h1 = histogram(RWT,m1:binstep:m2,'Facecolor','b');
            h2 = histogram(RKO,m1:binstep:m2,'Facecolor','r');
            Z = max([h1.Values,h2.Values]);
            plot(nanmedian(RWT)*[1 1],[0 Z],'Color','b','Linewidth',2);
            plot(nanmedian(RKO)*[1 1],[0 Z],'Color','r','Linewidth',2);
            plot_hist_signif(pvalue,nanmedian(RWT),nanmedian(RKO),Z,Z,0,Z/10);
            axis tight;
            xlabel('Correlations'); ylabel('Number of cells');
            title([titM{m},': ',titC{py}]);
            
            % PLOT MEAN PER ANIMAL ------------------------------------
            RWT = cellfun(@nanmedian, rwt);                                   % Average per animal
            RKO = cellfun(@nanmedian, rko);
            
            if py == 2                                                      % If computing mean over interneurons per animal
                RWT(8) = [];                                                % REMOVE 8th WT ANIMAL. IT HAS ONLY 1 INTERNEURON
                RWT(3) = [];                                                % REMOVE 3rd WT ANIMAL. IT HAS ONLY 2 INTERNEURONS
                RKO(1) = [];                                                % REMOVE 1st KO ANIMAL. IT HAS ONLY 2 INTERNEURONS
            end
            
            [pvalue, testtype] = significance(RWT,RKO,'unequal'); testtype;
            pv = ones(2,2);
            pv(1,2) = pvalue * 3;                                           % BONFERRONI CORRECTION OVER 2 CELL TYPE COMPARISONS
            
            subplot(2,9,(m-1)*9+(py-1)*3+3);hold on;
            plot_mean_SE(1:2,{RWT;RKO},pv,cols);
            set(gca,'Xtick',1:2,'XTickLabel',titKO);
            title('Mean per animal');
            xlim([0.4 2.6]);
        end
        
        
        % PY-IN CORRELATIONS
        rwt = cRho(KO==1,x,m);                                              % Keep WT mice
        rko = cRho(KO==2,x,m);                                              % And KO mice
        
        % POOL ALL CELLS ------------------------------------------
        RWT = cell2mat(rwt);                                                % Concatenate cells from all animals
        RKO = cell2mat(rko);
        
        RWT(isnan(RWT)) = [];
        RKO(isnan(RKO)) = [];
        
        [pvalue, testtype] = significance(RWT,RKO,'unequal'); testtype;
        pvalue = pvalue * 3;                                                % BONFERRONI CORRECTION OVER 2 CELL TYPE COMPARISONS
        
        subplot(2,9,(m-1)*9+(py-1)*6+[1 2]);hold on;
        m1 = corrlims(m,1);
        m2 = corrlims(m,2);
        binstep = (m2-m1)/nbins;
        h1 = histogram(RWT,m1:binstep:m2,'Facecolor','b');
        h2 = histogram(RKO,m1:binstep:m2,'Facecolor','r');
        Z = max([h1.Values,h2.Values]);
        plot(nanmedian(RWT)*[1 1],[0 Z],'Color','b','Linewidth',2);
        plot(nanmedian(RKO)*[1 1],[0 Z],'Color','r','Linewidth',2);
        plot_hist_signif(pvalue,nanmedian(RWT),nanmedian(RKO),Z,Z,0,Z/10);
        axis tight;
        xlabel('Correlations'); %ylabel('Number of cells');
        title([titM{m},': PY-IN']);
        
        % PLOT MEAN PER ANIMAL ------------------------------------
        RWT = cellfun(@nanmedian, rwt);                                       % Average per animal
        RKO = cellfun(@nanmedian, rko);
        
        RWT(8) = [];                                                        % REMOVE 8th WT ANIMAL. IT HAS ONLY 1 INTERNEURON
        RWT(3) = [];                                                        % REMOVE 3rd WT ANIMAL. IT HAS ONLY 2 INTERNEURONS
        RKO(1) = [];                                                        % REMOVE 1st KO ANIMAL. IT HAS ONLY 2 INTERNEURONS
                
        [pvalue, testtype] = significance(RWT,RKO,'unequal'); testtype;
        pv = ones(2,2);
        pv(1,2) = pvalue * 3;                                               % BONFERRONI CORRECTION OVER 3 CELL TYPE COMPARISONS
        
        subplot(2,9,(m-1)*9+(py-1)*6+3);hold on;
        plot_mean_SE(1:2,{RWT;RKO},pv,cols);
        set(gca,'Xtick',1:2,'XTickLabel',titKO);
        title('Mean per animal');
        xlim([0.4 2.6]);
    end
end


