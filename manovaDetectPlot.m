%%Brendan Thorn - 03/16
% manovaDetectPlot.m
%
% This code is essentially the same as anovaDetectPlot.m, but uses a
% multi-variate one-way ANOVA instead of the univariate version

fprintf('Beginning manovaDetectPlot.m at time:\n');
c=clock; fprintf('%.0f:%.0f\n',c(4),c(5));

% Provide folder name for the correct data
folder = 'M20100302_645';
% Load timing info and LFP data structure into workspace for trial info
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Set type of trials to inspect
type = 'sacs';
% Find number of, and particular, appropriate trials based on 'type'
% NOTE: THIS PART MUST SELECT SAME TRIALS AS anova_code4.m !!!
validTrials = 1:length(info);
if strcmp(type,'sacs')
    validSacs = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133 & validTrials'~=302);
    numSacs = length(validSacs); numPurs = 0;
    numTrainTrials = numSacs;
elseif strcmp(type,'purs')
    validPurs = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    numPurs = length(validPurs); numSacs = 0;
    numTrainTrials = numPurs;
elseif strcmp(type,'both')
    validSacs = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133 & validTrials'~=302);
    validPurs = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    numSacs = length(validSacs);
    numPurs = length(validPurs);
    % If we use both and want to train with same number of trials, use this
    numTrainTrials = min(numSacs,numPurs);
end
numTrials = max(numSacs,numPurs);

% ****** TO SET ****************************************************
% Set which channels to look at
chanType = 'mip';
if strcmp(chanType,'mip')
    valChannels = 1:36;
elseif strcmp(chanType,'pmd')
    valChannels = 37:48;
end
channels = 1:length(valChannels);
numChannels = length(channels);

% Set (approximate) bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15; 
bounds(3) = 25; bounds(4) = 35; 
bounds(5) = 45; bounds(6) = 55; 
bounds(7) = 65; bounds(8) = 75; 
bounds(9) = 85; bounds(10) = 95; 
bounds(11) = 105;bounds(12) = 125;
bounds(13) = 150;
numFreqBands = length(bounds);

% Load relevant channel spectrograms into workspace
% NOTE: spects matrix is cropped concatenation of spectrograms from all
% channels. The time dimension (3) begins at index 201 of the full
% spectrogram matrices, so to acces idx 281, use idx 81 in spect. (This
% note is only true for the original (mip and) pmdLookSpects.mat files).
if strcmp(chanType,'mip')
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    load(handle);
elseif strcmp(chanType,'pmd')
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    load(handle);
else % Load relevant channels (all for now)
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    mip=load(handle);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    t = mip.t; f = mip.f;
    clear mip pmd;
end

% TO SET
% Choose which times to use for imminent eye movement
moveIdxs = 78;
% Choose which times to use for NO imminent eye movement (It is okay to 
% have moveIdxs inside noMoveIdxs as it is illustrative of no difference)
noMoveIdxs = [50:100]';
numPeriods = size(noMoveIdxs,1);
times = t(50:100);

sacDetectPVals = zeros(numChannels,numPeriods);
purDetectPVals = zeros(numChannels,numPeriods);

% Perform anova for the eye movements
for chan = 1:numChannels
    % Get log-spectrogram info for correct channel
    spect = s(spects(channels(chan),:,:,:));

    % Average over desired frequency bands
    [freqBands, fBounds] = create_freqBands(spect, f, bounds);

    if ~strcmp(type,'purs')
        trainIdxs = randperm(numSacs,numTrainTrials);
        trainSacs = validSacs(trainIdxs);
        % Check for tuning over specified time steps (near eye movement onset)
        for l = 1:numPeriods
            bandVals = reshape(freqBands(:,[moveIdxs,noMoveIdxs(l,1)],trainSacs),numFreqBands,[])';
            moveVals = repmat([ones(length(moveIdxs),1);zeros(length(noMoveIdxs(l,:)),1)],numTrainTrials,1);
            [~,sacDetectPVals(chan,l)] = manova1(bandVals,moveVals);
        end
    end
    
    if ~strcmp(type,'sacs')
        % Purs next
        trainIdxs = randperm(numPurs,numTrainTrials);
        trainPurs = validPurs(trainIdxs);
        % Check for tuning over specified time steps (near eye movement onset)
        for l = 1:numPeriods
            bandVals = reshape(freqBands(:,[moveIdxs,noMoveIdxs(l,1)],trainPurs),numFreqBands,[])';
            moveVals = repmat([ones(length(moveIdxs),1);zeros(length(noMoveIdxs(l,:)),1)],numTrainTrials,1);
            [~,purDetectPVals(chan,l)] = manova1(bandVals,moveVals);
        end
    end

    fprintf('Completed LOO ANOVA for channel %d.\n',valChannels(chan));
end

% Set p-value that determines if a channel is significantly tuned
pVal = 0.05;

sacSigPVals = s(sum(sacDetectPVals<pVal,1));
purSigPVals = s(sum(purDetectPVals<pVal,1));

handle = sprintf('../analysis/workspaces/%s/%sDetMANOVA_plotData_140316.mat',folder,chanType);
%save(handle,'sacDetectPVals','purDetectPVals','sacSigPVals','purSigPVals','times','folder','chanType','fBounds','pVal','validSacs','validPurs','moveIdxs','noMoveIdxs');

%%%%%%%%%%%% Plot of number of temporally-tuned channels vs time %%%%%%%%%%%%%%

figure;
plot(times(1:end),sacSigPVals);

% Insert vertical line at last time used before movement onset and at
% movement onset itself
line([t(moveIdxs(1)) t(moveIdxs(1))], [-.5 36],'LineWidth',1);line([t(moveIdxs(end)) t(moveIdxs(end))], [-.5 36],'LineWidth',1); line([0 0], [-.5 36],'LineWidth',1); 
line([t(moveIdxs(1)-1) t(moveIdxs(1)-1)], [-.5 36],'LineWidth',1,'Color','r');line([t(moveIdxs(end)+1) t(moveIdxs(end)+1)], [-.5 36],'LineWidth',1,'Color','r');
xlim([times(1) times(end)]);ylim([-0.1 36]);
title('MIP ANOVA - Temporal Tuning for Saccade Detection by Frequency Band (p < 0.05)');
xlabel('Time Relative to Saccade [s]');ylabel('Number of Significantly Tuned Channels /36');
