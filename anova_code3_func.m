function [] = anova_code3_func( validTrials, type, pVal, confLevel, folder )
%%Brendan Thorn - 07/15
% anova_code3_func.m
%
% This is the same program as anova_code2.m, but accepts input args so it
% can be called from other programs
%
% Runs (one-way) anova on several frequency bands of LFP spectrograms
% to determine if directional tuning exists.
%
% This is a refactored version of anova_code.m.

fprintf('Beginning anova_code3_func.m!\n');

handle = sprintf('../data/%s/TT.mat',folder);
load(handle);

% Load an LFP data structure into workspace for trial info
filename = sprintf('../data/%s/LFPChan01.mat',folder);
load(filename);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

group = info(:,4);
dir = info(:,1)-128;

% Find trials within valid trials for each type of eye movement
purs = validTrials(group(validTrials)==66 & ~isnan(TT.T66on(validTrials)) & dir(validTrials)>0 & dir(validTrials)<5);
sacs = validTrials(group(validTrials)==65 & ~isnan(TT.T65sac(validTrials)) & dir(validTrials)>0 & dir(validTrials)<5);
fixs = validTrials(group(validTrials)==67 & ~isnan(TT.T8(validTrials)) & dir(validTrials)>0 & dir(validTrials)<5);
% saccades which are farther than median from reach
sacThresh = median((TT.T8(sacs)-TT.T65sac(sacs)));
isoSacs = (TT.T8-TT.T65sac) > sacThresh;

% ****** TO SET ****************************************************
% Set trials to whichever eye movement type you are interested in *****
trials1 = sacs;
trials2 = purs;

% Set which channels to look at
channels = 1:48;
numChannels = length(channels);
% Set which time steps to look at
timeIdxs = 271:2:309;
numTimeSteps = length(timeIdxs);
% Set bounds for frequency bands to look at
bounds(1) = 6; bounds(2) = 12; 
bounds(3) = 18; bounds(4) = 24; 
bounds(5) = 30; bounds(6) = 36; 
bounds(7) = 42; bounds(8) = 48; 
bounds(9) = 54; bounds(10) = 66; 
bounds(11) = 80;bounds(12) = 100;
bounds(13) = 150;
numFreqBands = length(bounds);

sacTuningPVals = zeros(numChannels,numFreqBands,numTimeSteps);
purTuningPVals = zeros(numChannels,numFreqBands,numTimeSteps);


%for each channel
for chan = 1:numChannels
    % TO SET
    % Get spectrogram/coherogram/cross-power for channel
    handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,channels(chan),47);
    [spect,t,f] = getCohs(1,handle,channels(chan),channels(1),folder);
    % Center time vector
    t=t-7.5;
    fprintf('Spectrogram retrieved for channel %d.\n',channels(chan));
    
%     % (OPTIONAL) Normalize the spectrogram
%     % Normalize spectrogram to activity ~800-500ms before eye movement
%     % First find correct timing
%     numNormSteps = 3;
%     %idxTime = (TT.T6-alignTimes)/1000;
%     idxTime = repmat(-.675,numTrials,1);
%     idxs = zeros(numTrials,numNormSteps);
%     for i = 1:numTrials
%         idxs(i,1) = find((idxTime(i,:) < t),1);
%         idxs(i,2:end) = idxs(i,1)+1:idxs(i,1)+numNormSteps-1;
%     end
    
    % Normalize the spectrogram
    %normSpect = normalize_spect(log(permute(spect,[2 1 3])),idxs);
    % Or don't normalize, but pretend you did so code still works
    normSpect = log(permute(spect,[2 1 3]));
    % Get rid of infs created by log to avoid problems later
    normSpect(isinf(normSpect)) = 0;
    
    % Average over desired frequency bands
    % Separate normalized spectrograms into bands between 0-200 Hz. Avoid noisy areas ie.
    % ~60 Hz and harmonics
    [freqBands, fBounds] = create_freqBands(normSpect, f, bounds);

    % Check for tuning over specified time steps (near eye movement onset)
    for k = 1:numFreqBands
        for j = 1:numTimeSteps
            if strcmp(type,'sacs')
                sacTuningPVals(chan,k,j) = anovan(s(freqBands(k,timeIdxs(j),trials1)),dir(trials1,1));
                close(1);
            elseif strcmp(type,'purs')
                purTuningPVals(chan,k,j) = anovan(s(freqBands(k,timeIdxs(j),trials2)),dir(trials2,1));
                close(1);
            elseif strcmp(type,'both')
                sacTuningPVals(chan,k,j) = anovan(s(freqBands(k,timeIdxs(j),trials1)),dir(trials1,1));
                close(1);
                purTuningPVals(chan,k,j) = anovan(s(freqBands(k,timeIdxs(j),trials2)),dir(trials2,1));
                close(1);
            end
        end
    end
    times=t(timeIdxs);
    fprintf('Done channel %d of %d!\n',chan,numChannels);
end

% Create some dummy values in case variables aren't used later (because
% they will still be saved)
sacSigPVals = 0; purSigPVals = 0;
sacSigChans = 0; purSigChans = 0;
if strcmp(type,'sacs')
    sacSigPVals = s(sum(sacTuningPVals(1:numChannels,:,:)<pVal,1));
    sacSigChans = cell(numFreqBands,numTimeSteps);
    for i = 1:numFreqBands
        for j = 1:numTimeSteps
            sacSigChans{i,j}=find(sacTuningPVals(1:numChannels,i,j)<pVal);
        end
    end
elseif strcmp(type,'purs')
    purSigPVals = s(sum(sacTuningPVals(1:numChannels,:,:)<pVal,1));
    purSigChans = cell(numFreqBands,numTimeSteps);
    for i = 1:numFreqBands
        for j = 1:numTimeSteps
            purSigChans{i,j}=find(purTuningPVals(1:numChannels,i,j)<pVal);
        end
    end
elseif strcmp(type,'both')
    sacSigPVals = s(sum(sacTuningPVals(1:numChannels,:,:)<pVal,1));
    sacSigChans = cell(numFreqBands,numTimeSteps);
    purSigPVals = s(sum(purTuningPVals(1:numChannels,:,:)<pVal,1));
    purSigChans = cell(numFreqBands,numTimeSteps);
    for i = 1:numFreqBands
        for j = 1:numTimeSteps
            sacSigChans{i,j}=find(sacTuningPVals(1:numChannels,i,j)<pVal);
            purSigChans{i,j}=find(purTuningPVals(1:numChannels,i,j)<pVal);d
        end
    end
end

%handle = sprintf('../analysis/workspaces/%s/%sTunedChans_%d_NN_jan1816.mat',folder,type,confLevel);
%save(handle,'sacTuningPVals','purTuningPVals','sacSigChans','purSigChans','sacSigPVals','purSigPVals','times','timeIdxs','fBounds');

fprintf('Completed anova_code3_func.m!\n');
end