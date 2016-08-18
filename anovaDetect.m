%%Brendan Thorn - 02/16
% anovaDetect.m
%
% This code aims to determine if certain frequebcy bands of LFP
% spectrograms carry a significant amount of information about eye movement
% timing. Initially, it will be modelled after anova_code4.m, which is used
% to determine which frequency bands carry tuning information. This code
% should be compatible with svm_LookDetect2.m.

fprintf('Beginning anovaDetect.m at time:\n');
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
type = 'both';
% Find number of, and particular, appropriate trials based on 'type'
% NOTE: THIS PART MUST SELECT SAME TRIALS AS anova_code4.m !!!
validTrials = 1:length(info);
if strcmp(type,'sacs')
    validSacs = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    numSacs = length(validSacs); numPurs = 0;
elseif strcmp(type,'purs')
    validPurs = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    numPurs = length(validPurs); numSacs = 0;
elseif strcmp(type,'both')
    validSacs = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    validPurs = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    numSacs = length(validSacs);
    numPurs = length(validPurs);
end
numTrials = max(numSacs,numPurs);

% ****** TO SET ****************************************************
% Set which channels to look at
chanType = 'both';
channels = 1:48;
numChannels = length(channels);
% Set a default p-value (this isn't important if you save lots of vars)
pVal = 0.01;
% Set (approximate) bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15;
bounds(3) = 25; bounds(4) = 35;
bounds(5) = 45; bounds(6) = 55;
bounds(7) = 65; bounds(8) = 75;
bounds(9) = 85; bounds(10) = 95;
bounds(11) = 105;bounds(12) = 125;
bounds(13) = 150;
numFreqBands = length(bounds);

sacDetectPVals = zeros(numChannels,numSacs,numFreqBands);
purDetectPVals = zeros(numChannels,numPurs,numFreqBands);

% Load relevant channel spectrograms into workspace
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

% Choose which times to use for imminent eye movementt
moveIdxs = 78:79;
% Create vector that classifies each time step as "MOVE"=1 or "NO MOVE"=2
move = zeros(length(t),1);
move(moveIdxs) = 1;
% Choose which times to use for ANOVA
tIdxs = [41:2:63 64:75 moveIdxs];

% Perform anova for the eye movements
for chan = 1:numChannels
    % Get log-spectrogram info for correct channel
    spect = s(spects(channels(chan),:,:,:));

    % Average over desired frequency bands
    [freqBands, fBounds] = create_freqBands(spect, f, bounds);

    % Now loop through trials. When trial = i, we leave out the i-th trial
    % so it can be used for cross-val later. Sacs first
    for trial = 1:numSacs
        trainSacs = [validSacs(1:trial-1) validSacs(trial+1:end)];
        numTrainTrials = length(trainSacs);
        % Check for tuning over specified time steps (near eye movement onset)
        for k = 1:numFreqBands
                bandVals = reshape(s(freqBands(k,tIdxs,trainSacs)),[],1);
                moveVals = repmat(move(tIdxs),numTrainTrials,1);
                sacDetectPVals(chan,trial,k) = anovan(bandVals,moveVals,'display','off');
        end
    end
    % Purs next
    for trial = 1:numPurs
        trainPurs = [validPurs(1:trial-1) validPurs(trial+1:end)];
        numTrainTrials = length(trainPurs);
        % Check for tuning over specified time steps (near eye movement onset)
        for k = 1:numFreqBands
                bandVals = reshape(s(freqBands(k,tIdxs,trainPurs)),[],1);
                moveVals = repmat(move(tIdxs),numTrainTrials,1);
                purDetectPVals(chan,trial,k) = anovan(bandVals,moveVals,'display','off');
        end
    end
    fprintf('Completed LOO ANOVA for channel %d.\n',channels(chan));
end
% At this point, we have the 3D sacDetectPVals matrix, which tells us for 
% each channel, and for each trial left out for each channel, the p-value 
% at which the frequency bands are temporally tuned. Use this matrix to 
% make to create two 2D matrices that 1) give the number of tuned channels 
% for each frequency band, for each left out trial; and 2) list exactly 
% which channels are tuned and comprise the number given in the previous 
% matrix. Compare to arbitrarily chosen pVal.
sacSigPVals = s(sum(sacDetectPVals<pVal,1));
sacSigChans = cell(numSacs,numFreqBands);
for trial = 1:numSacs
    for fb = 1:numFreqBands
        sacSigChans{trial,fb}=find(sacDetectPVals(:,trial,fb)<pVal);
    end
end

purSigPVals = s(sum(purDetectPVals<pVal,1));
purSigChans = cell(numPurs,numFreqBands);
for trial = 1:numPurs
    for fb = 1:numFreqBands
        purSigChans{trial,fb}=find(purDetectPVals(:,trial,fb)<pVal);
    end
end

times=t(tIdxs);

handle = sprintf('../analysis/workspaces/%s/LOO_detANOVA_NN_170316_7879.mat',folder);
save(handle,'sacDetectPVals','purDetectPVals','sacSigChans','purSigChans','sacSigPVals','purSigPVals','times','tIdxs','fBounds','pVal','validSacs','validPurs','moveIdxs');

