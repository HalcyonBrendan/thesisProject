%%Brendan Thorn - 03/16
% anovaDetectPlot.m
%
% This code borrows from anovaDetect.m to produce a plot that demonstrates
% whether band-times at various times throughout the trial are
% significantly different from the band-times very close to the movement.

fprintf('Beginning anovaDetectPlot.m at time:\n');
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
    valChannels = 1:32;
elseif strcmp(chanType,'pmd')
    valChannels = 33:48;
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
    handle = sprintf('../data/%s/mipLookSpects_240_024.mat',folder);
    load(handle);
elseif strcmp(chanType,'pmd')
    handle = sprintf('../data/%s/pmdLookSpects_240_024.mat',folder);
    load(handle);
else % Load relevant channels (all for now)
    handle = sprintf('../data/%s/mipLookSpects_240_024.mat',folder);
    mip=load(handle);
    handle = sprintf('../data/%s/pmdLookSpects_240_024.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    t = mip.t; f = mip.f;
    clear mip pmd;
end

% TO SET
% Choose which times to use for imminent eye movement
moveIdxs = [61 63];
% Create an nxm matrix of time indices to compare to moveIdxs. Each row
% should represent a "period" to compare with ANOVA to the "move period"
% defined by moveIdxs
noMoveIdxs = [21:73]';
numPeriods = size(noMoveIdxs,1);
times = t(21:73);

sacDetectPVals = zeros(numChannels,numFreqBands,numPeriods);
purDetectPVals = zeros(numChannels,numFreqBands,numPeriods);

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
        for k = 1:numFreqBands
            for l = 1:numPeriods
                bandVals = reshape(s(freqBands(k,[moveIdxs,noMoveIdxs(l,:)],trainSacs)),[],1);
                moveVals = repmat([ones(length(moveIdxs),1);zeros(length(noMoveIdxs(l,:)),1)],numTrainTrials,1);
                sacDetectPVals(chan,k,l) = anovan(bandVals,moveVals,'display','off');
            end
        end
    end
    
    if ~strcmp(type,'sacs')
        % Purs next
        trainIdxs = randperm(numPurs,numTrainTrials);
        trainPurs = validPurs(trainIdxs);
        % Check for tuning over specified time steps (near eye movement onset)
        for k = 1:numFreqBands
            for l = 1:numPeriods
                bandVals = reshape(s(freqBands(k,[moveIdxs,noMoveIdxs(l,:)],trainPurs)),[],1);
                moveVals = repmat([ones(length(moveIdxs),1);zeros(length(noMoveIdxs(l,:)),1)],numTrainTrials,1);
                purDetectPVals(chan,k,l) = anovan(bandVals,moveVals,'display','off');
            end
        end
    end

    fprintf('Completed LOO ANOVA for channel %d.\n',valChannels(chan));
end

% Set p-value that determines if a channel is significantly tuned
pVal = 0.01;

sacSigPVals = s(sum(sacDetectPVals<pVal,1));
purSigPVals = s(sum(purDetectPVals<pVal,1));

handle = sprintf('../analysis/workspaces/%s/%sDetANOVA_plotData_240_024_170316.mat',folder,chanType);
%save(handle,'sacDetectPVals','purDetectPVals','sacSigPVals','purSigPVals','times','folder','chanType','fBounds','pVal','validSacs','validPurs','moveIdxs','noMoveIdxs');

%%%%%%%%%%%% Plot of number of temporally-tuned channels vs time %%%%%%%%%%%%%%

% Get some colours to use for lines in plot
colors = rand(13,3);
%=[0.0539481928176814,0.548419485731466,0.949370243451058;0.880415938951547,0.0126686341217107,0.647537226254477;0.614411014426690,0.987721046933028,0.0382624597507391;0.00342308456076834,0.0374395360441805,0.123355202520449;0.963402881241132,0.194147564660471,0.293973863482836;0.644568825739417,0.440747029351203,0.680998914337276;0.488862294275121,0.243484263440880,0.644282796413461;0.178710699560875,0.662681633681037,0.187649735135656;0.388925444913291,0.519870901392384,0.416803601643293;0.610385431090712,0.720587737805813,0.782348041459805;0.652039122948022,0.104973057264752,0.112291674642342;0.853291109133284,0.200642858627722,0.413806386178684;0.480308365434127,0.946990212518994,0.776160594864731]

bands = 1:13; 
figure;
for i = 1:length(bands)
    plot(times(1:end), sacSigPVals(bands(i),:),'color',colors(i,:));
    hold on;
end
hold off;
% Insert vertical line at last time used before movement onset and at
% movement onset itself
line([-.2 -.2], [-.5 32],'LineWidth',1);line([-.1 -.1], [-.5 32],'LineWidth',1); line([0 0], [-.5 32],'LineWidth',1); 
xlim([times(1) times(end)]);ylim([-0.1 32]);
title('MIP ANOVA - Temporal Tuning for Saccade Detection by Frequency Band (p < 0.01)');
xlabel('Time Relative to Saccade [s]');ylabel('Number of Significantly Tuned Channels /32');
