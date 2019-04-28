function [S, lout] = remove_bad_units(animal,S)

%% SET THE BAD UNITS
if strcmp(animal,'M2420')
    badsets = [2];
    badunits = {4};
    
elseif strcmp(animal,'M2422')
    badsets = [4];
    badunits = {6};
    
elseif strcmp(animal,'M2453')
    badsets = [1 2];
    badunits = {[3]; [1]};
    
elseif strcmp(animal,'M2455')
    badsets = [2 3 6];
    badunits = {[1]; [1 2 5]; [6]};
    
elseif strcmp(animal,'M2730')
    badsets = [2];
    badunits = {[3]};
    
elseif strcmp(animal,'M2870')
    badsets = [1 3 6 7];
    badunits = {[6]; [5]; [7]; [3]};
    
elseif strcmp(animal,'MD1A')
    badsets = [1 3];
    badunits = {[1]; [7]};
    
elseif strcmp(animal,'MD1B')
    badsets = [2 3 5];
    badunits = {[6]; [4]; [5]};
    
elseif strcmp(animal,'MD2A')
    badsets = [4 5];
    badunits = {[1]; [7]};
    
elseif strcmp(animal,'MD3A')
    badsets = [1 2];
    badunits = {[1 7]; [6]};
    
elseif strcmp(animal,'MD4B')
    badsets = [3];
    badunits = {[2]};
    
else
    badsets = {};
    badunits = {};
end

lout = cellfun(@length, badunits);
lout = sum(lout);

%% FIND THE ACTUAL INDICES OF THE BAD SETS
sudirfiles = subdir(fullfile(['../',animal],'*.spike.mat'));  % Make list with all .spike.mat files
ls = length(sudirfiles);
for st = 1:ls                % For each Set file
    sufile = sudirfiles(st).name;           % Get its name
    setnum(st) = str2double(sufile(strfind(sufile,'Set')+3)); % STORE FILE NAME WITH SET INFO
end
for s = 1:length(badsets)
    badsets(s) = find(setnum == badsets(s));
end

%% REMOVE BAD UNITS
for i = 1:length(badsets)
    S.Peak{badsets(i)}(badunits{i}) = [];
    S.Trough{badsets(i)}(badunits{i}) = [];
    S.PTdist{badsets(i)}(badunits{i}) = [];
    S.PTamp{badsets(i)}(badunits{i}) = [];
    S.Halfwidth{badsets(i)}(badunits{i}) = [];
    
    S.Waveforms{badsets(i)}(badunits{i},:) = [];
    
    S.mRate{badsets(i)}(:,:,badunits{i}) = [];
    S.BIndex{badsets(i)}(:,:,badunits{i}) = [];
    S.CSIndex{badsets(i)}(:,:,badunits{i}) = [];
    
    S.Phases{badsets(i)}(:,:,:,badunits{i}) = [];
    S.mPhases{badsets(i)}(:,:,:,badunits{i}) = [];
    S.PhaseVec{badsets(i)}(:,:,:,badunits{i}) = [];
    S.PhaseRt{badsets(i)}(:,:,:,badunits{i}) = [];
end
