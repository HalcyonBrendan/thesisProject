%%Brendan Thorn - 10/15
% lfp_tuningCurves.m
% 
% Builds (4-direction) tuning "curves" for pursuit and saccade eye
% movements. The curves are actually just mean power for movement
% in each direction for given freq band at given time in trials.
%
% Note: Before running, load LFP data (eg. LFPChan01.mat - only need
% one channel to get trial info) and strobe timing data (eg. TT.mat)
% into workspace. Also adjust handles within code as necessary.

%some variables and structures we will need later
numChans = 48;
numDirs = 4;

memTimes = -70:4:90;
memTimes(2,:) = memTimes(1,:)*.025;
numTimes = length(memTimes);

group = Events.dataSummary(:,4);
direction = Events.dataSummary(:,1)-128;

purs = group==66;
fixs = group==67;
sacs = group==65;
% saccades which are farther than median from reach
sacThresh = median((TT.T8(sacs)-TT.T65sac(sacs)));
isoSacs = (TT.T8-TT.T65sac) > sacThresh;

% ****** TO BE SET ****************************************************
% set trials to whichever eye movement type you are interested in *****
trials = purs;

% declare arrays that contain mean power info for each direction - one for
% each frequency band
tunCurves0010 = zeros(numChans, numTimes, numDirs);
tunCurves1020 = zeros(numChans, numTimes, numDirs);
tunCurves2030 = zeros(numChans, numTimes, numDirs);
tunCurves3040 = zeros(numChans, numTimes, numDirs);
tunCurves4060 = zeros(numChans, numTimes, numDirs);

%for each channel
for i = 1:numChans
    
    %every 4th channel, load next four channels
    if mod(i,4)==1
        %handle = sprintf('../analysis/workspaces/M20100301_246/chanArray%d-%dAug0315.mat',i,i+3);
        handle = sprintf('../analysis/workspaces/monkey_data/chanArray%d-%dJul1715.mat',i,i+3);
    end
    clear chanArray;
    load(handle);
    
    for j = 1:numTimes
        %consider mean power in ~0-10Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(1:5,221-6+memTimes(1,j),trials)));
        %now get mean power for trials of interest
        for k = 1:numDirs
            tunCurves0010(i,j,k) = mean(temp(direction(trials)==k));
        end
        %consider mean power in ~10-20Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(6:11,221-6+memTimes(1,j),trials)));
        for k = 1:numDirs
            tunCurves1020(i,j,k) = mean(temp(direction(trials)==k));
        end
        %consider mean power in ~20-30Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(12:16,221-6+memTimes(1,j),trials)));
        for k = 1:numDirs
            tunCurves2030(i,j,k) = mean(temp(direction(trials)==k));
        end 
        %consider mean power in ~30-40Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(17:21,221-6+memTimes(1,j),trials)));
        for k = 1:numDirs
            tunCurves3040(i,j,k) = mean(temp(direction(trials)==k));
        end
        %consider mean power in ~40-60Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(22:30,221-6+memTimes(1,j),trials)));
        for k = 1:numDirs
            tunCurves4060(i,j,k) = mean(temp(direction(trials)==k));
        end
    end
    disp(i);
    
end

% create tuning curve array curvesxxyyzz for memTimes[xx], for freq band
% yy-zz Hz
curves224060 = normr(squeeze(tunCurves4060(:,22,:)));
for i =1:48
    rowSum=sum(curves224060(i,:));
    for j=1:4
        curves224060(i,j)=curves224060(i,j)/rowSum;
    end
end

curves332030 = normr(squeeze(tunCurves2030(:,33,:)));
for i =1:48
    rowSum=sum(curves332030(i,:));
    for j=1:4
        curves332030(i,j)=curves332030(i,j)/rowSum;
    end
end

curves252030 = normr(squeeze(tunCurves2030(:,25,:)));
for i =1:48
    rowSum=sum(curves252030(i,:));
    for j=1:4
        curves252030(i,j)=curves252030(i,j)/rowSum;
    end
end
