function Get_Motion(animal)

%% SET DIRECTORIES
outputdir = fullfile('..','Analysis Results',animal);
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end
motfile = fullfile(outputdir,'Motion.mat');

%% SET PARAMETERS
Fs = 25e3;
Fsd = 1e3;

initime = 5*60;                                                             % Iinital time is after the first 5 min

critV1 = 0.8;                                                               % Criterion for movement detection (ratio of velocity mean)
critV2 = 5e-3;                                                              % Criterion of immobility
if strcmp(animal(1:2),'MD')
    critV1 = 2;   
    critV2 = 0.016;
end

%% READ OXYTOCIN POINT
[oxypoint,animals] = xlsread('Oxt_datapoints.xlsx',1,'','basic');

k = contains(animals,animal);
oxypoint = oxypoint(k);
oxytime = oxypoint/Fs;

%% DETECT MOVEMENT
Vfile = fullfile(['../',animal],'stimuli','stimuli_run.mat');
load(Vfile);
V = stimuli.stim;
V = V(1 : Fs/Fsd : end);
clear stimuli
[Move,MoveT,V,Vtime,critV] = detect_movement(V,Fsd,[critV1,critV2]);

%% COMPUTE PRE-POST OXYTOCIN MOTION AND IMMOBILITY DURATION
sudirfiles = subdir(fullfile(['../',animal],'*.spike.mat'));                % Make list with all .spike.mat files
L = load(sudirfiles(1).name);                                               % Load FIRST SET
spiketimes = double(L.spike_time)/1e6;                                      % spike times in sec

z = L.cluster_id .* spiketimes;
z = z(:);
lastspike = max(z);                                                         % APPROXIMATE TOTAL RECORDING TIME (NOT USING LFP DURATION SINCE IT MAY CONTAIN NOISE IN THE END, e.g. M2422)

if strcmp(animal,'MD2B') 
    lastspike = Vtime(end);                                                 % EXCEPTION becauser 'lastspike' is very early
end

MoveT = MoveT(MoveT(:,1) > initime,:);                                      % Exclude first 5 mins
MoveT = MoveT(MoveT(:,2) <= lastspike,:);                                   % Exclude time after last unit spike

Mpre = find(MoveT(:,2) < oxytime - 5*60,1,'last');                          % Find last motion segment ENDING pre-Oxytocin
Mpost = find(MoveT(:,1) > oxytime + 5*60,1,'first');                        % Find first motion segment STARTING post-Oxytocin
MoveTox{1} = MoveT(1:Mpre,:);
MoveTox{2} = MoveT(Mpost:end,:);

time_Ox(1) = oxytime - 5*60;                                                % Time duration preOxytocin
time_Ox(2) = lastspike - (oxytime + 5*60);                                  % Time duration postOxytocin
for x = 1:2
    total_time_mot(x,1) = sum(diff(MoveTox{x},[],2));                       % Total Motion-time (pre or post Oxytocin)
    total_time_mot(x,2) = time_Ox(x) - total_time_mot(x);                   % Total Immobility-time (pre or post Oxytocin)
end
disp(' ');
disp(['Total time of motion before oxytocin = ', num2str(total_time_mot(1,1)),' sec']);
disp(['Total time of immobility before oxytocin = ', num2str(total_time_mot(1,2)),' sec']);
disp(' ');
disp(['Total time of motion after oxytocin = ', num2str(total_time_mot(2,1)),' sec']);
disp(['Total time of immobility after oxytocin = ', num2str(total_time_mot(2,2)),' sec']);
disp(' ');

%% PLOT MOTION
screensize = get(0,'Screensize');
fm = figure('Name',['Motion profile of animal ',animal]);
set(gcf, 'PaperUnits', 'points', 'Units', 'points');
set(gcf, 'Position', [10 0 screensize(3)*0.9 screensize(4)*0.7]);

subplot(1,1,1);hold on;
plot(Vtime(~Move), V(~Move),'-k','Linewidth',1.5);
plot(Vtime(Move), V(Move),'-b','Linewidth',1.5);
line(oxytime*[1 1],[0 max(V)],'Color','k','LineWidth',1.5);
line((oxytime-5*60)*[1 1],[0 max(V)],'Color','k','LineStyle','--','LineWidth',1.5);
line((oxytime+5*60)*[1 1],[0 max(V)],'Color','k','LineStyle','--','LineWidth',1.5);
line(initime*[1 1],[0 max(V)],'Color','g','LineWidth',1.5);
line(lastspike*[1 1],[0 max(V)],'Color','r','LineWidth',1.5);
line([Vtime(1) Vtime(end)],critV(1)*[1 1],'Color','k','LineWidth',1.5);
line([Vtime(1) Vtime(end)],critV(2)*[1 1],'Color','k','LineWidth',1.5);
xlabel('time(sec)'); ylabel('Motion'); xlim([Vtime(1),Vtime(end)]);ylim([0 max(V)]);box on;
drawnow;

%% SAVE MOTION DATA AND FIGURE
save(motfile,'V','MoveT','total_time_mot','MoveTox','lastspike');

saveas(fm,fullfile(outputdir,'Motion.fig'));
saveas(fm,fullfile(outputdir,'Motion.jpg'));

