function LFP_bandpower(animal,channel,Frange)

%% SET DIRECTORIES AND LOAD MOTION DATA
outputdir = fullfile('..','Analysis Results',animal);
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end
outputfile = fullfile(outputdir,['LFPpower_ch',num2str(channel),'.mat']);

LFP1000dir = fullfile(['../',animal],'LFP','LFP1000');
S = load(fullfile('..','Analysis Results',animal,'Motion.mat'));     % Load animal motion file

%% SET PARAMETERS AND FILTERS
Fs = 25e3;
Fsd = 1e3;
lfr = size(Frange,1);

%% READ OXYTOCIN POINT
[oxypoint,animals] = xlsread('Oxt_datapoints.xlsx',1,'','basic');
k = contains(animals,animal);
oxypoint = oxypoint(k);                                                     % Find the correct oxy datapoint
oxytime = oxypoint/Fs;                                                      % turn to timepoint (using the original sampling rate)

%% COMPUTE PRE-POST OXYTOCIN MOTION AND IMMOBILITY LFP-INDEXES
dt = 1/Fsd;
time = 0 : dt : S.lastspike;                                                % Time vector for the whole recording
lt = length(time);

indexes = cell(2,2);
for x = 1:2                                                                 % For each oxy condition
    indexes{x,1} = zeros(1,lt);                                             % Zeros vector for each time point                
    M = S.MoveTox{x};                                                       % Keep corresponding motion segments
    for j = 1:size(M,1)                                                     % For each motion segment
        [~,k1] = min(abs(time - M(j,1)));                                   % Find closest time points to segment edges
        [~,k2] = min(abs(time - M(j,2)));
        indexes{x,1}(k1:k2) = 1;                                            % Turn zeros to ones for segments of motion
    end
    indexes{x,2} = ~indexes{x,1};                                           % Keep opposite segments for immobility OVER THE ENTIRE SESSION
end

k = find(time <= 5*60,1,'last');                                            % Find last index before first 5 min pass
indexes{1,2}(1:k) = 0;                                                      % Set pre-oxy immobility indexes up to that index back to zero
k = find(time <= oxytime-5*60,1,'last');                                    % Find last index before oxytime (minus 5 min)
indexes{1,2}(k+1:end) = 0;                                                  % Set pre-oxy immobility indexes beyond that index back to zero

k = find(time >= oxytime+5*60,1,'first');                                   % Find first index after oxytime (plus 5 min)
indexes{2,2}(1:k) = 0;                                                      % Set post-oxy immobility indexes up to that index back to zero
k = find(time <= S.lastspike,1,'last');                                     % Find last index before last spike
indexes{2,2}(k+1:end) = 0;                                                  % Set post-oxy immobility indexes beyond that index back to zero

%% POWER SPECTRUM OF LFP FROM GIVEN CHANNEL
disp(['Filtering LFP from channel ' num2str(channel)]);
LFPfile = fullfile(LFP1000dir,['LFPvoltage_ch' num2str(channel) '.mat']);
L = load(LFPfile);
LFPvoltage = double(L.LFPvoltage);
clear L

LFPspec = cell(2,2);
sg = 10;

for x = 1:2                                                                 % For each oxy condition
    for m = 1:2                                                             % For each motion condition
        LFPind = LFPvoltage(logical(indexes{x,m}));                         % Keep the LFP at the corrseponding segments
        if ~isempty(LFPind)
            [sp,~,~] = sp2a2_m1(0,LFPind',LFPind',Fsd,sg);                  % Compute power spectrum
            LFPspec{x,m} = sp(:,2) + log10(2*pi);                            
            
            if any(isinf(LFPspec{x,m}))
                LFPspec{x,m} = [];
            end
        end
    end
end
freq = sp(:,1);

%% COMPUTE POWER
LFPpower = zeros(2,2,lfr);
for x = 1:2                                                                 % For each oxy condition
    for m = 1:2                                                             % For each motion condition
        for f = 1:lfr                                                       % For each frequency range
            LFPind = LFPvoltage(logical(indexes{x,m}));                     % Keep the LFP during corresponding segments
            if ~isempty(LFPind)
                LFPpower(x,m,f) = bandpower(LFPind,Fsd,[Frange(f,1), Frange(f,2)]); % Compute power
            end
        end
    end
end

%% SAVE
save(outputfile,'LFPpower','LFPspec','indexes','freq');

