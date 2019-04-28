function Motion_Analysis(animals,KO)

%% MOTION META-ANALYSIS
la = length(animals);

Vm = cell(la,1);
Vi = cell(la,1);
Vf = cell(la,1);
mV = cell(la,1);
totV = cell(la,1);
Tm = cell(la,1);
totTm = zeros(la,2);

Fs = 1e3;
fTm = zeros(1:la);

for a = 1:la
    S = load(fullfile('..','Analysis Results',animals{a},'Motion.mat'));    % Load animal unit file
    x = 1;                                                                  % Analyze only pre-oxy segments
    
    Tm{a} = S.MoveTox{x};
    totTm(a,:) = 100 * S.total_time_mot(x,:) / sum(S.total_time_mot(x,:));  % Percentage of time in motion/immobility over total recording time
    fTm(a) = size(Tm{a},1) / sum(S.total_time_mot(x,:));                    % Frequency of motion bouts (#bouts / recording duration)
    
    Vm{a} = cell(length(Tm{a}),1);          
    t = 0 : 1/Fs : (length(S.V)-1)/Fs;                                      % Time axis over total recording
    for i = 1:length(Tm{a})                                                 % For each motion segment
        Vm{a}{i} = S.V(t >= Tm{a}(i,1) &  t <= Tm{a}(i,2));                 % Keep velocity during it
        %         Vm{a}{i} = Vm{a}{i} / max(Vm{a}{i});
        
        if Tm{a}(i,1)-1 > 0 &  Tm{a}(i,1)+1 < t(end)                        % If there is a full second before and after the motion initiation
            [~,k] = min(abs(t - Tm{a}(i,1)));                               % Find the timepoint closest to motion initiation 
            Vi{a}(i,:) = S.V(k-1000 : k+1000)';%/S.V(k+1000);               % Keep velocity around it (+- 1 sec)
        end
        if Tm{a}(i,2)-1 > 0 &  Tm{a}(i,2)+1 < t(end)                        % Same for motion termination
            [~,k] = min(abs(t - Tm{a}(i,2)));
            Vf{a}(i,:) = S.V(k-1000 : k+1000)';%/S.V(k-1000)';
        end
    end
    
    mV{a} = cellfun(@mean,Vm{a});                                           % Keep mean velocity PER SEGMENT
    totV{a} = cellfun(@trapz, Vm{a});% ./ diff(Tm{a},[],2);                 % Keep velocity integral (distance) PER SEGMENT
end
clear S

%% PLOT
pv = ones(2,2);

figure;
subplot(241);
[pv(1,2),ttype] = significance(fTm(KO == 1),fTm(KO == 2),'unequal');ttype;
plot_mean_SE(1:2,{fTm(KO == 1),fTm(KO == 2)},pv,['b','r']);
ylabel('Mean frequency of motion bouts per animal (#bouts/recording time)'); 

subplot(242);
dTm = cellfun(@(x) diff(x,[],2), Tm,'UniformOutput',0);                     % Get motion bout durations
dTm = cellfun(@mean,dTm);                                                   % Get mean per animal
pv = ones(2,2);
pv(1,2) = significance(dTm(KO == 1),dTm(KO == 2),'unequal');
plot_mean_SE(1:2,{dTm(KO == 1),dTm(KO == 2)},pv,['b','r']);
ylabel('Mean duration of motion bouts per animal (sec)'); 

subplot(243);
pv = ones(2,2);
[pv(1,2),ttype] = significance(totTm(KO==1,1),totTm(KO==2,1),'unequal');ttype;
plot_mean_SE(1:2,{totTm(KO==1,1),totTm(KO==2,1)},pv,['b','r'])
ylabel('% time in motion per animal'); 

subplot(244);
pv = ones(2,2);
[pv(1,2),ttype]  = significance(totTm(KO==1,2),totTm(KO==2,2),'unequal');ttype;
plot_mean_SE(1:2,{totTm(KO==1,2),totTm(KO==2,2)},pv,['b','r'])
ylabel('% time in immobility per animal'); 

subplot(245);
mv = cellfun(@mean,mV);                                                     % Get mean average-velocity per animal
pv = ones(2,2);
[pv(1,2),ttype] = significance(mv(KO == 1),mv(KO == 2),'unequal');ttype;
plot_mean_SE(1:2,{mv(KO == 1),mv(KO == 2)},pv,['b','r']);
ylabel('Mean velocity per motion bout per animal (a.u.)'); 

subplot(246);
totv = cellfun(@mean,totV);
pv = ones(2,2);
[pv(1,2),ttype] = significance(totv(KO == 1),totv(KO == 2),'unequal');ttype;
plot_mean_SE(1:2,{totv(KO == 1),totv(KO == 2)},pv,['b','r']);
ylabel('Mean distance per motion bout per animal (a.u.)'); 

% Plot motion initiation
subplot(247);
hold on;
mViWT = cellfun(@(x) mean(x,1), Vi(KO == 1),'UniformOutput',0);
mViWT = cell2mat(mViWT);
mViKO = cellfun(@(x) mean(x,1), Vi(KO == 2),'UniformOutput',0);
mViKO = cell2mat(mViKO);
for a = 1:sum(KO == 1)
    plot(-1:0.001:1,mViWT(a,:),'b')
end
for a = 1:sum(KO == 2)
    plot(-1:0.001:1,mViKO(a,:),'r')
end
fill_plot(-1:0.001:1,mean(mViWT,1),0,std(mViWT,[],1)/sqrt(sum(KO == 1)),'b')
plot(-1:0.001:1,mean(mViWT,1),'b','Linewidth',2);
fill_plot(-1:0.001:1,mean(mViKO,1),0,std(mViKO,[],1)/sqrt(sum(KO == 2)),'r')
plot(-1:0.001:1,mean(mViKO),'r','Linewidth',2);
xlabel('time (sec)'); ylabel('Mean velocity during motion initiation per animal(a.u.)')

pv = ones(1,2001);
ttype = zeros(1,2001);
for i = 1:2001
    [pv(i),ttype(i)] = significance(mViWT(:,i),mViKO(:,i),'unequal');
end
ttype;
[~,pv] = fdr(pv);
plot(find(pv < 0.05),max(mean(mViWT,1))*ones(sum(pv < 0.05),1),'*k');

% Plot motion termination
subplot(248);
hold on;
mVfWT = cellfun(@(x) mean(x,1), Vf(KO == 1),'UniformOutput',0);
mVfWT = cell2mat(mVfWT);
mVfKO = cellfun(@(x) mean(x,1), Vf(KO == 2),'UniformOutput',0);
mVfKO = cell2mat(mVfKO);
for a = 1:sum(KO == 1)
    plot(-1:0.001:1,mVfWT(a,:),'b')
end
for a = 1:sum(KO == 2)
    plot(-1:0.001:1,mVfKO(a,:),'r')
end
fill_plot(-1:0.001:1,mean(mVfWT,1),0,std(mVfWT,[],1)/sqrt(sum(KO == 1)),'b')
plot(-1:0.001:1,mean(mVfWT,1),'b','Linewidth',2);
fill_plot(-1:0.001:1,mean(mVfKO,1),0,std(mVfKO,[],1)/sqrt(sum(KO == 2)),'r')
plot(-1:0.001:1,mean(mVfKO),'r','Linewidth',2);
xlabel('time (sec)'); ylabel('Mean velocity during motion termination per animal (a.u.)')   

pv = ones(1,2001);
ttype = zeros(1,2001);
for i = 1:2001
    [pv(i),ttype(i)] = significance(mVfWT(:,i),mVfKO(:,i),'unequal');
end
ttype;
[~,pv] = fdr(pv);
plot(find(pv < 0.05),max(mean(mVfWT,1))*ones(sum(pv < 0.05),1),'*k');
