function Get_Rates(animals,binlen,spikethr)

%% GET PARAMETERS
load(fullfile('..','Analysis Results','Clusters.mat'));

la = length(animals);

initime = 5*60;                                                             % Iinital time is after the first 5 min

% Read oxytocin point
[oxy,anim] = xlsread('Oxt_datapoints.xlsx',1,'','basic');
for a = 1:la                                                                % For each animal
    k = contains(anim,animals(a));
    oxypoint = oxy(k);
    oxytime = oxypoint/1000;
    ftime(a) = oxytime - 5*60;
end

Rate = cell(la,2,2,2);
mRate = cell(la,2,2,2);
sRate = cell(la,2,2,2);
ISI = cell(la,2,2,2);
totspikes = cell(la,2,2,2);

%% COMPUTE RATES
for a = 1:la                                                                % For each animal
    a
    S = load(fullfile('..','Analysis Results',animals{a},'Units.mat'));     % Load animal unit file
    load(fullfile('..','Analysis Results',animals{a},'Motion.mat'),'MoveTox'); % Load animal unit file
    
    S = remove_bad_units(animals{a},S);                                     % REMOVE PRE-SELECTED BAD UNITS
    
    ls = length(S.mRate);                                                   % Number of sets for that animal
    
    for x = 1%:2                                                            % For each condition
        Mt = MoveTox{x};                                                    % Keep motion segements there
        
        for py = 1:2                                                        % For each cell type
            count = 1;                                                      % Start cell counter
            for st = 1:ls                                                   % For each set
                cells = (clustered_units{a}{st} == py);                     % Get indexes of that type of cells
                lc = sum(cells);
                
                ST = S.SpikeTimes{st}(cells);                               % Keep spiketimes of only those cells
                
                for u = 1:lc                                                % For each cell
                    % MOTION RATES ----------------------------------------
                    STu = ST{u}{x,1};                                       % Keep only spikes during motion
                    Ru = [];
                    isi = [];
                    for m = 1:size(Mt,1)                                    % For each motion segment
                        R = rate_per_segment(STu,Mt(m,1),Mt(m,2),binlen);   % Compute firing rates within the segment
                        Ru = cat(1,Ru,R');                                  % Concatenate with rates of previous motion segments
                        z = STu(STu >= Mt(m,1) & STu <= Mt(m,2));
                        isi = cat(1,isi,diff(z));
                    end
                    Rate{a,py,x,1}(:,count) = Ru;                           % Store cell's mean rate over all motion segments
                    mRate{a,py,x,1}(count,1) = mean(Ru);                         % Store cell's mean rate over all motion segments
                    sRate{a,py,x,1}(count,1) = var(Ru);                         % Store cell's variance over all motion segments
                    ISI{a,py,x,1}(count,1) = mean(isi);
                    totspikes{a,py,x,1}(count) = length(STu);
                    
                    % IMMOBILITY RATES ------------------------------------
                    STu = ST{u}{x,2};                                   % Keep only spikes during immobility
                    isi = [];
                    if ~ isempty(STu)
                        Ru = rate_per_segment(STu,initime,Mt(1,1),binlen); % Compute firig rates from first immobility spike to first motion segment
                        R = rate_per_segment(STu,Mt(end,2),ftime(a),binlen); % Compute firig rates from end of last motion segment to last immobility spike
                        Ru = cat(1,Ru',R');                               % Concatenate with rates of previous motion segments
                        
                        z = STu(STu >= initime & STu <= Mt(1,1));
                        isi = cat(1,isi,diff(z));
                        z = STu(STu >= Mt(end,2) & STu <= ftime(a));
                        isi = cat(1,isi,diff(z));
                        
                        for m = 1:size(Mt,1)-1                                % For each motion segment (except last)
                            R = rate_per_segment(STu,Mt(m,2),Mt(m+1,2),binlen); % Compute firig rates between the segment and the next one
                            Ru = cat(1,Ru,R');                               % Concatenate with rates of previous motion segments
                            
                            z = STu(STu >= Mt(m,1) & STu <= Mt(m,2));
                            isi = cat(1,isi,diff(z));
                        end
                        Rate{a,py,x,2}(:,count) = Ru;
                        mRate{a,py,x,2}(count,1) = mean(Ru);                         % Store cell's mean rate over all motion segments
                        sRate{a,py,x,2}(count,1) = var(Ru);                         % Store cell's variance over all motion segments
                        ISI{a,py,x,2}(count,1) = mean(isi);
                        totspikes{a,py,x,2}(count) = length(STu);
                    end
                    count = count + 1;
                end
            end
            Rate{a,py,x,1} = sparse(Rate{a,py,x,1});
            Rate{a,py,x,2} = sparse(Rate{a,py,x,2});
        end
    end
end

%% COMPUTE CORRELATIONS
Rho = cell(la,2,2,2);
cRho = cell(la,2,2);

for a = 1:la                                                                % For each animal
    a
    for x = 1%:2
        for m = 1:2
            for py = 1:2
                R = Rate{a,py,x,m};                                         % Keep its rates
                Rho{a,py,x,m} = corr(R,R,'Type','Pearson');                 % Compute correlations between all pairs of cells of same type
                
                k = find(totspikes{a,py,x,m} < spikethr);
                Rho{a,py,x,m}(k,:) = nan;
                Rho{a,py,x,m}(:,k) = nan;
            end
            cRho{a,x,m} = corr(Rate{a,1,x,m},Rate{a,2,x,m},'Type','Pearson'); % Compute correlations between all pairs of PY-IN
            
            k = find(totspikes{a,1,x,m} < spikethr);
            cRho{a,x,m}(k,:) = nan;
            k = find(totspikes{a,2,x,m} < spikethr);
            cRho{a,x,m}(:,k) = nan; 
        end
    end
end

%% SAVE
save(fullfile('..','Analysis Results',['Rates_',num2str(binlen),'.mat']),'Rate','mRate','sRate','ISI','Rho','cRho','totspikes');
