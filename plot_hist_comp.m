function plot_hist_comp(A,KO,fname,cols,nbins)
 
titOx = {'Pre-Oxy', 'Post-Oxy'};
titM = {'Motion', 'Immobility'};
titC = {'PY','IN'};
titKO = {'WT','KO'};

if nargin == 4
    nbins = 30;                                                                 % Number of histogram bins
end

for x = 1%:2
    figure('Name',fname);

    for m = 1:2                                                             % For motion/immobility
        for py = 1:2                                                        % For each cell type
            rwt = A(KO==1,py,x,m);                                          % Keep WT mice
            rko = A(KO==2,py,x,m);                                          % And KO mice
            
            if size(rwt{8},2) == 0 % Fix one empty problematic correlation matrix
                rwt{8} = rwt{8}';
            end
            
            % POOL ALL CELLS ------------------------------------------
            RWT = cell2mat(rwt);                                            % Concatenate cells from all animals
            RKO = cell2mat(rko);
            
            RWT(isnan(RWT)) = [];
            RKO(isnan(RKO)) = [];
            
            [pvalue, testtype] = significance(RWT,RKO,'unequal'); testtype;
            pvalue = pvalue * 2;                                            % BONFERRONI CORRECTION OVER 2 CELL TYPE COMPARISONS
            
            subplot(2,6,(m-1)*6+(py-1)*3+[1 2]);hold on;
            m1 = min([RWT;RKO]);
            m2 = max([RWT;RKO]); 
%             m2=50  % For problematic ISI
            binstep = (m2-m1)/nbins;                                       
            h1 = histogram(RWT,m1:binstep:m2,'Facecolor','b');
            h2 = histogram(RKO,m1:binstep:m2,'Facecolor','r');
            Z = max([h1.Values,h2.Values]);
            plot(nanmean(RWT)*[1 1],[0 Z],'Color','b','Linewidth',2);
            plot(nanmean(RKO)*[1 1],[0 Z],'Color','r','Linewidth',2);
            plot_hist_signif(pvalue,nanmean(RWT),nanmean(RKO),Z,Z,0,Z/10);
            axis tight;
            xlabel(fname); ylabel('Number of cells');
            title([titOx{x},' ',titM{m},': ',titC{py}]);
            
            % PLOT MEAN PER ANIMAL ------------------------------------
            RWT = cellfun(@nanmean, rwt);                                   % Average per animal
            RKO = cellfun(@nanmean, rko);
%             
%             if py == 2                                                      % If computing mean over interneurons per animal
%                 RWT(8) = [];                                                % REMOVE 8th WT ANIMAL. IT HAS ONLY 1 INTERNEURON
%                 RWT(3) = [];                                                % REMOVE 3rd WT ANIMAL. IT HAS ONLY 2 INTERNEURONS
%                 RKO(1) = [];                                                % REMOVE 1st KO ANIMAL. IT HAS ONLY 2 INTERNEURONS
%             end
%             
            [pvalue, testtype] = significance(RWT,RKO,'unequal'); testtype;
            pv = ones(2,2);
            pv(1,2) = pvalue * 2;                                           % BONFERRONI CORRECTION OVER 2 CELL TYPE COMPARISONS
            
            subplot(2,6,(m-1)*6+(py-1)*3+3);hold on;
            plot_mean_SE(1:2,{RWT;RKO},pv,cols);
            set(gca,'Xtick',1:2,'XTickLabel',titKO);
            ylabel('Mean per mouse');
        end
    end
end