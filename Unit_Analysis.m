function Unit_Analysis(animal)

addpath('CircStat2012a');
screensize = get(0,'Screensize');

cols = [9 117 220;      %blue
        255 64 0]/255;	% red
    
%% SET DIRECTORIES
outputdir = fullfile('..','Analysis Results',animal);
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end
unitsfile = fullfile(outputdir,'Units.mat');

BSdir = fullfile(['../',animal],'LFP','BackSub');
LFP1000dir = fullfile(['../',animal],'LFP','LFP1000');

%% LOAD MOTION DATA
motfile = fullfile(outputdir,'Motion.mat');
load(motfile);

%% SET PARAMETERS AND FILTERS
Fs = 25e3;
Fsd = 1e3;

initime = 5*60;                     % Iinital time is after the first 5 min

Fi = 600;
Ff = 6000;
Wn = [Fi Ff]/(Fs/2); %low freq, high freq, sr/2
[Bb,Ba] = butter(1,Wn,'bandpass'); %finds the coefficients for the butter filter

window = 35;
time = (-window:window)/Fs * 1000;   %time axis in ms

bnd = {'Delta';'Theta';'Beta';'Sgam';'Fgam'};

Frange = [1,4;
    5,11;
    12,30;
    30,55;
    80,110];
lfr = size(Frange,1);
for f = 1:lfr
    Hd{f} = bandpass(Fs,Frange(f,1),Frange(f,2));
end

lr = 3+lfr;
lc = 7;

burstDT = 0.020;  % 20 ms

titOx = {'Pre-Oxy', 'Post-Oxy'};
titM = {'Motion', 'Immobility'};

%% READ OXYTOCIN POINT
[oxypoint,animals] = xlsread('Oxt_datapoints.xlsx',1,'','basic');
k = contains(animals,animal);
oxypoint = oxypoint(k);
oxytime = oxypoint/Fs;

%% GET CHANNEL SETS
load('probe128A.mat'); %contains probelayout
shank = 2;
sets = {1:7;             % electdoe sets of Shank2
    8:18;
    19:27;
    28:36;
    37:47;
    48:57;
    58:64};
for st = 1:7          % For each set of electrode LOCATIONS
    channels{st} = probelayout(sets{st},shank); % Keep the electrode IDs
end

%% SET UP MEASURES TO COMPUTE
sudirfiles = subdir(fullfile(['../',animal],'*.spike.mat'));  % Make list with all .spike.mat files
ls = length(sudirfiles);

Peak = cell(ls,1);
Trough = cell(ls,1);
PTdist = cell(ls,1);
PTamp = cell(ls,1);
Halfwidth = cell(ls,1);
Waveforms = cell(ls,1);
mRate = cell(ls,1);
BIndex = cell(ls,1);
CSIndex = cell(ls,1);
InterSpike = cell(ls,1);
Phases = cell(ls,1);
mPhases = cell(ls,1);
PhaseVec = cell(ls,1);
PhaseRt = cell(ls,1);
SpikeTimes = cell(ls,1);

%% LOAD MAT FILE DATA AND GET SPIKE WAVEFORMS
for st = 1:ls                                                               % For each Set file
    sufile = sudirfiles(st).name;                                           % Get its name
    setnum = sufile(strfind(sufile,'Set')+3);                               % and its number
    disp(['Working on waveforms from set ',num2str(setnum)]);
    maxSp = load(sufile);                                                   % Load it!
    
    clust = maxSp.cluster_id;                                               % for that spike file/set, stores each unit as a column
    spiketimes = double(maxSp.spike_time)/1e6;                              % spike times in sec of all units
    setnum = str2double(setnum);                                            % STORE FILE NAME WITH SET INFO
    chanset = channels{setnum};
    Nch = length(chanset);                                                  % Numbe of channels in the set
    
    Nsp = sum(clust,1);                                                     % Number of spikes per unit.
    clust(:,Nsp == 0) = [];                                                 % Remove clusters that have no spikes (!!)
    Nsp(Nsp == 0) = [];                                                 
    Nu = size(clust,2);                                                     % Final number of units
        
    %% GET SPIKETIMES OF EACH UNIT
    unitsp = cell(Nu,1);
    for u = 1:Nu                                                            %for each unit (which is a column and contains only 0 and 1s)
        unitsp{u} = spiketimes(clust(:,u) == 1);                            % Keep the corresponding spike times
    end
    
    spikeW = cell(Nu,Nch);
    Spikes = cell(Nu,1);
    Waveforms{st} = zeros(Nu,2*window+1);
    
    %% ISOLATE PERI-SPIKE LFP SEGMENTS
    parfor j = 1:Nch                                                        % for each channel in that set
        LFPfile = fullfile(BSdir,['LFPvoltage_ch' num2str(chanset(j)) '.mat']); % Get the corresponding LFP filename
        if exist(LFPfile,'file')
            disp(['Channel ' num2str(chanset(j)),' (',num2str(j) '/' num2str(length(chanset)),')'])
            L = load(LFPfile);                                              % Load the file
            for u = 1:Nu                                                    % for each unit
                spikeW{u,j} = zeros(Nsp(u),2*window + 1);                   % Allocate space for waveforms              
                startsamp = round(unitsp{u}*Fs);                            % compute the datapoint over the recording  of that spike
                for s = 1:Nsp(u)                                            %for every spike time
                    spikeW{u,j}(s,:) = L.LFPvoltage(startsamp(s) + (-window:window)); % Store waveform from lfp, starting 30 datapoints before and ending 30 after startstamp, which is when unit is identified?
                end
            end
        else
            disp(['Channel ' num2str(chanset(j)),' (',num2str(j) '/' num2str(length(chanset)),') BAD CHANNEL, SKIPPING']); % If it doesnt exist skip
        end
    end
    clear L
    
    
    %% GET CHANNEL WITH LARGEST WAVEFORM AND KEEP
    mspike = cellfun(@mean, spikeW,'UniformOutput',0);                      % Mean waveform of each unit
    mspike = cellfun(@abs, mspike,'UniformOutput',0);                       % The absolute value of the mean waveform
    spikeamps = cellfun(@max, mspike);                                      % Peak of waveform
    [~,bestchans] = max(spikeamps,[],2);                                    % Find set channel where that is maximum
    
    Peak{st} = zeros(Nu,1);
    Trough{st} = zeros(Nu,1);
    PTdist{st} = zeros(Nu,1);
    PTamp{st} = zeros(Nu,1);
    Halfwidth{st} = zeros(Nu,1);
    Waveforms{st} = zeros(Nu,2*window+1);
    
    K = zeros(Nu,3);
    Z = zeros(Nu,2);
    
    for u = 1:Nu                                                            % For each unit in the set
        maxSp = spikeW{u,bestchans(u)};                                     % Keep the waveforms on channel with maximum waveform peak
        for s = 1:Nsp(u)                                                    % For each spike of the unit
            maxSp(s,:) = filtfilt(Bb,Ba, maxSp(s,:));                       % Filter waveform (sampled at 25kHz) at 600-6000 Hz
        end
        Spikes{u} = -maxSp;                                                 % Reverse filtered waveforms
        W = mean(-maxSp,1);                                                 % and keap mean waveform
        
        [A,k] = findpeaks(abs(W),'SortStr','descend');                      % Get peaks in descending order and their times
        [K(u,:),order] = sort(k(1:3));                                      % Sort the times of the largest 3 peaks
        A = A(1:3);                                                         % Keep peaks
        A = A(order).*sign(W(K(u,:)));                                      % Order them by time and switch to their original sign
        
        Peak{st}(u) = A(2);                                                 % Waveform peak (second peak)
        Trough{st}(u) = A(3);                                               % Waveform trough (third peak)
        PTdist{st}(u) = (K(u,3)-K(u,2))/Fs * 1000;                          % Peak - Trough distance in ms
        PTamp{st}(u) = A(2) - A(3);                                         % Peak-Trough amplitude
        
        [~,z] = findpeaks(-abs(W - Peak{st}(u)/2),'SortStr','descend');     % Find times closest to the half-peak points (essentially in their original order)
        Z(u,:) = z(1:2);                                                    % Keep the two closest
        Halfwidth{st}(u) = abs(Z(u,2) - Z(u,1))/Fs * 1000;                  % Halfwidth in ms
        
        Waveforms{st}(u,:) = W;                                             % Store waveforms
    end
    
    
    %% SPLIT SPIKES
    SpikeTimes{st} = cell(Nu,1);
    
    for u = 1:Nu                                                            % For each unit in the set
        SpikeTimes{st}{u} = cell(2,2);
        S = unitsp{u};                                                      % Keep its spiketimes
        for x = 1:2                                                         % For each oxy-condition
            if x == 1                                                       % If pre-oxy
                Sox = S(S > initime & S < oxytime - 5*60);                  % Keep spikes after initial time and before 5min-preOxy
            else
                Sox = S(S > oxytime + 5*60 & S <= lastspike);               % Or spikes after 5min-postOxy and before 'lastspike'
            end
            
            k = [];
            for m = 1:size(MoveTox{x},1)                                    % For each motion segment
                sm = find(Sox >= MoveTox{x}(m,1) & Sox < MoveTox{x}(m,2));  % Find spikes within it
                k = [k; sm];                                            
            end
            SpikeTimes{st}{u}{x,1} = Sox(k);                                    % Keep motion spikes
            k = setdiff(1:length(Sox), k);                                  % Get all the remaining spikes (immobility)
            SpikeTimes{st}{u}{x,2} = Sox(k);                                    % Store immobility spikes
        end
    end
    
    
    %% MEAN FIRING RATE, BURST INDEX, CSI
    mRate{st} = zeros(2,2,Nu);
    BIndex{st} = zeros(2,2,Nu);
    CSIndex{st} = zeros(2,2,Nu);
    InterSpike{st} = zeros(2,2,Nu);
    
    for u = 1:Nu                                                            % For each unit
        for x = 1:2                                                         % And each condition
            for m = 1:2
                S = SpikeTimes{st}{u}{x,m};                                 % Keep its spiketimes
                ns = length(S);                                             % Number of spikes
                
                mRate{st}(x,m,u) = ns / total_time_mot(x,m);                % Mean rate
                
                Dtsp = diff(S);                                             % Spike time distance
                InterSpike{st}(x,m,u) = mean(Dtsp);                         % Mean inter-spike interval
                
                bursts = [0; Dtsp < burstDT];                               % Get bursty spikes
                BIndex{st}(x,m,u) = sum(bursts) / ns * 100;                 % Burst index = percentage of spikes closer than 20 ms

                amp = zeros(1,ns);
                allamps = max(Spikes{u},[],2);                              % Get the spike amplitudes of all spikes
                for s = 1:ns
                    amp(s) = allamps(unitsp{u} == S(s));                    % Keep amplitudes of corresponding spikes only
                end               
                bi = find(bursts - circshift(bursts,-1) == -1);             % Find zeros followed by a 1 (1st spike in burst)
                bf = find(bursts - circshift(bursts,-1) == 1);              % Find ones followed by a 0 (last spike in burst)
                decr = zeros(1,length(bi));
                incr = zeros(1,length(bi));
                for b = 1:length(bi)                                        % For each burst
                    bamps = amp(bi(b) : bf(b));                             % Keep spike amplitudes in the burst
                    dbamps = diff(bamps);                                   % Compute change in amplitude relative to predecessor (excludes 1st spike in burst)
                    decr(b) = sum(dbamps < 0);                              % Count amplitude drops
                    incr(b) = sum(dbamps > 0);                              % And increases
                end
                CSIndex{st}(x,m,u) = (sum(decr) - sum(incr)) / sum(bf-bi) * 100; % Ratio of total number of amplitude drops - increases over total number of bursty spikes
            end
        end
    end
    
    %% PLOT WAVERFORMS
    for u = 1:Nu
        figs{u} = figure('Name',['Animal ', animal,'. Set ',num2str(setnum),'. Unit ',num2str(u)]);
        set(gcf, 'PaperUnits', 'points', 'Units', 'points');
        set(gcf, 'Position', [10 0 screensize(3)*0.9 screensize(4)*0.7]);
        
        subplot(lr,lc,repmat((0:lc:(2+lfr)*lc)',1,3)+(1:3)); hold on;
        for s = 1:100:Nsp(u)                                                % Plot every 20th spike
            [~,maxp] = max(Spikes{u}(s, K(u,2) +(-5:5)));                   % Find the spike's peak point close to the mean peak point
            maxp = 5-maxp+1;                                             % Find their distance
            maxp = maxp/Fs*1000;                                            % Turn to time distance
            plot(time+maxp,Spikes{u}(s,:),'Color',[0.5 0.5 0.5]);
        end
        plot(time,Waveforms{st}(u,:),'k','Linewidth',2);
        line([time(1),time(end)],[0 0],'Color','k');
        xlim([time(11) time(61)]);                                          % Show only -1:1 sec
        title(['LFP channel ',num2str(chanset(bestchans(u)))]);
        xlabel('time (ms)');
        ylabel(['Unit ',num2str(u)]);
        
        plot(time(K(u,2)),Peak{st}(u),'^r','Linewidth',1.5);
        plot(time(K(u,3)),Trough{st}(u),'^r','Linewidth',1.5);
        line(time(K(u,2))*[1 1], [0 Peak{st}(u)],'Color','r','Linewidth',1.5);
        line(time(K(u,3))*[1 1], [0 Trough{st}(u)],'Color','r','Linewidth',1.5);
        line(time([K(u,2) K(u,3)]), Peak{st}(u)*[1.2 1.2],'Color','k','Linewidth',1.5);
        plot(time([Z(u,1) Z(u,2)]), Peak{st}(u)/2*[1 1],'Color','g','Linewidth',1.5);
        box on;
        
        ax{1} = subplot(lr,lc,[4 11]);
        h = bar3(mRate{st}(:,:,u));
        h(1).FaceColor = 'b';h(2).FaceColor = 'r';
        zlabel('Mean Rate (Hz)');box on;axis tight;
        
        ax{2} = subplot(lr,lc,[5 12]);
        h = bar3(BIndex{st}(:,:,u));
        h(1).FaceColor = 'b';h(2).FaceColor = 'r';
        zlabel('Burst Index');box on;axis tight;
        
        ax{3} = subplot(lr,lc,[6 13]);
        h = bar3(CSIndex{st}(:,:,u));
        h(1).FaceColor = 'b';h(2).FaceColor = 'r';
        zlabel('CS Index');box on;axis tight;
        for x = 1:3
            ax{x}.XTickLabel = titM;                                        % (counter-intuitive but correct)
            ax{x}.YTickLabel = {'Pre','Post'};
        end
    end
    drawnow;
    
    %% GET PHASES
    disp(['Analyzing spike phases from set ',num2str(setnum)]);
    Phases{st} = cell(lfr,2,2,Nu);
    mPhases{st} = zeros(lfr,2,2,Nu);
    PhaseVec{st} = zeros(lfr,2,2,Nu);
    PhaseRt{st} = ones(lfr,2,2,Nu);
    
    for u = 1:Nu                                                            % For each unit in the set
        disp(['Filtering LFP from channel ' num2str(chanset(bestchans(u)))]);
        LFPfile = fullfile(LFP1000dir,['LFPvoltage_ch' num2str(chanset(bestchans(u))) '.mat']);
        L = load(LFPfile);                                                  % Load the LFP that gave largest waveform (subsampled to 1kHz)
        L = double(L.LFPvoltage);
        
        for f = 1:lfr                                                       % For each frequency range
            LFPbp = filtfilt(Hd{f}.sosMatrix,Hd{f}.ScaleValues,L);          % Filter the LFP 
            y = hilbert(LFPbp);                                             % Compute Hilbert Transform 
            phase = angle(y);                                               % Get phase in radians (oscillation peaks corresponds to 0)
            for x = 1:2                                                     % For each oxy-condition
                for m = 1:2                                                 % For each motion-condition
                    Phases{st}{f,x,m,u} = phase(round(SpikeTimes{st}{u}{x,m}*Fsd))'; % Keep phases of correponding spike times (as column)
                    Phases{st}{f,x,m,u} = wrapToPi(Phases{st}{f,x,m,u}+pi); % Wrap in -pi - pi cycle
                    
                    if ~isempty(Phases{st}{f,x,m,u})
                        mPhases{st}(f,x,m,u) = circ_mean(Phases{st}{f,x,m,u}); % Mean phase
                        PhaseVec{st}(f,x,m,u) = circ_r(Phases{st}{f,x,m,u});   % Length of mean vector
                        PhaseRt{st}(f,x,m,u) = circ_rtest(Phases{st}{f,x,m,u});% Rayleigh test for non-uniformity
                    end
                end
            end
        end
    end
%     PhaseRt{st} = PhaseRt{st}*numel(PhaseRt{st});       % BONFERRONI CORRECTION!!
    
    %% PLOT PHASES
    a = linspace(-180,525);
    for u = 1:Nu
        figure(figs{u});
        for f = 1:lfr
            for x = 1:2
                subplot(lr,lc,24+(f-1)*lc+(x-1)*2+[1:2]); hold on;
                Mh = zeros(1,2);
                for m = 1:2
                    PH = radtodeg(Phases{st}{f,x,m,u});
                    PH = [PH; PH+360];                                      % Repeat all phases after adding a full cycle to show 2 oscillation cycles
                    if ~isempty(PH)
%                         histogram(PH,-180:5:540,'EdgeColor','none','FaceColor',cols(m,:));
                        K = histcounts(PH,-180:15:540);
                        K = K / sum(K) * 100;                               % Turn to percentage of spikes
                        bar(-180:15:525,K,'EdgeColor','none','FaceColor',cols(m,:),'FaceAlpha',0.5);
%                         Mh(m) = max(histc(PH,-180:5:540));
                        Mh(m) = max(K);
                    end
                end
                for m = 1:2
                    h = line(radtodeg(mPhases{st}(f,x,m,u))*[1 1], [0 1.2*max(Mh)],'Linewidth',1.5);
                    if m == 1
                        h.Color = 'b';
                    else
                        h.Color = 'r';
                    end
                    if PhaseRt{st}(f,x,m,u) * lfr * 2 < 0.05                % BONFERRONI CORRECTION OVER PLOT'S COMPARISONS
                        h.LineStyle = '-';                                  % FOR EACH UNIT SEPARATELY (NOT INCLUDING PRE-POST OXY CORRECTION)
                    else
                        h.LineStyle = ':';
                    end
                end
                plot(a, cos(deg2rad(a)+pi)*max(Mh)/2 + max(Mh)/2,'k','Linewidth',1.5);
                
                if f == lfr, xlabel('Phases');end;
                if f == 1, title(titOx{x}); end
                if x == 1, ylabel([bnd{f},'(',num2str(Frange(f,1)),'-',num2str(Frange(f,2)),'Hz)']);end
                set(gca,'Xtick',[-180:90:525],'Xticklabel',{'180';'270';'0';'90';'180';'270';'0';'90'});
                axis tight; box on;
            end
        end
        drawnow;
    end
    
    %% SAVE FIGURES AND CLOSE THEM
%     for u = 1:Nu
%         saveas(figs{u},fullfile('..','Analysis Results',animal,['Set ',num2str(setnum),' Unit ',num2str(u)]),'fig');
%         saveas(figs{u},fullfile('..','Analysis Results',animal,['Set ',num2str(setnum),' Unit ',num2str(u)]),'jpg');
%         close(figs{u});
%     end
end

save(unitsfile,'Peak','Trough','PTdist','PTamp','Halfwidth','Waveforms',...
    'SpikeTimes','mRate','BIndex','CSIndex','InterSpike',...
    'Phases','mPhases','PhaseVec','PhaseRt');

