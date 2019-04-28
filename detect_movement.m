function [Move,move_t,V,time,critV] = detect_movement(V,Fs,critV)

mindur = 0.5;                                                               % Minimum accepted duration of a segment
mindist = 0.5;                                                              % Minimum accepted time distance between two segments

N = length(V);
time = (0:N-1)/Fs;                                                          % Timescale for LFP (ALWAYS STARTS FROM ZERO. POSITION TIME STARTS LATER)

Wn = 1/(Fs/2);                                                              % lowpas at 1Hz
[b,a] = butter(1,Wn,'low');                                                 % finds the coefficients for the butter filter

if size(V,1) == 1
    V = V';
end

%% CALCULATE VELOCITY
V = V - mode(V);
V = abs(V);
V = filtfilt(b,a,V);                                                        % Lowpass the velocity
V(V<0) = 0;
V(1) = 0;
V(end) = 0;
critV(1) = critV(1)*mean(V);                                                    % Movement detection criterion

%% DETECT ALL SEGMENTS ABOVE THE THRESHOLD AND THEIR LIMITS
cross = (V >= critV(1));                                                      % Find the array locations over the threshold
mi = find(cross - circshift(cross,-1) == -1); % Find zeros preceding a 1
mi = mi + 1;
mf = find(cross - circshift(cross,-1) == 1);  % Find ones followed by a 0
rN = size(mi,1);

%% DETECT EDGES OF MOTION SEGMENTS
movei = zeros(rN,1);
movef = zeros(rN,1);
for i = 1:rN                                                                % For each segment
    movei(i) = find(V(1 : mi(i)) < critV(2),1,'last');                        % Keep the last location before segment-beginnning below the limit threshold
    movef(i) = find(V(mf(i) : end) < critV(2),1,'first') + mf(i) - 1;          % and the first location after segment-end below the limit threshold
end
move = [movei movef];

%% KEEP UNIQUE SEGMENTS
[~,m] = unique(move(:,1), 'rows');                                          % Keep unique events by checking repeating segment beginning entries
move = move(m,:);
[~,m] = unique(move(:,2), 'rows');                                          % Keep unique events by checking repeating segment end entries
move = move(m,:);

move_t = time(move);                                                        % Get time points from array locations

disp('Initially :')
[~,rdur] = motion_stats(move_t);                                             % Report

%% REMOVE SHORT AND OVERLAPPING SEGMENTS
rshort = (rdur < mindur);                                                   % If any segment lasts less than 0.1 sec discard it
move_t(rshort,:) = [];
move(rshort,:) = [];                                                        % Remove corresponding entries as well

k = move_t(2:end,1) - move_t(1:(end-1),2);                                  % Find the time distances between all segments
k = (k < mindist);                                                          % Find when the distance is not at least mindist (covers overlaps too)
k1 = [k; 0];                                                                % Keep locations of the preceding segment
k2 = [0; k];                                                                % and of the following segment
move(k1==1,2) = move(k2==1,2);                                              %  Unite them (use the ending of the following segment as end of the preceding one...
move(k2==1,:) = [];
move_t(k1==1,2) = move_t(k2==1,2);                                              %  Unite them (use the ending of the following segment as end of the preceding one...
move_t(k2==1,:) = [];                                                         %and delete the times of the following one)

disp('After short and close-segments removal: ')
rN = motion_stats(move_t);                                                   % Report

Move = zeros(1,length(V));
for t = 1:rN
    Move(move(t,1):move(t,2)) = 1;
end
Move = logical(Move);
%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------


function [rN,rdur] = motion_stats(mot)

rN = size(mot,1);                                                           % Number of motion segments in this session
rdur = mot(:,2)-mot(:,1);                                                   % and durations
rmax = max(rdur);                                                           % Maximum duration
disp([num2str(rN),' running sessions of average duration ', num2str(mean(rdur)*1000),...
    ' msec and maximum duration ',num2str(rmax*1000),' msec',' (SD = ',num2str(std(rdur)*1000),')'])

