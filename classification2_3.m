%%Brendan Thorn - 01/16
% classification2_3.m
%
% 25/01/16 - This is a refactored version of classification2_2.m that is
% designed to implement leave-one-out cross validation, as outlined in
% Scherberger et al. 2005. Cortical ...
%
% This is a reorganized version of classification2.m that requires
% less reading from disk to reduce runtime. The ordering (and indexing) is
% much less intuitive, but is necessary for practicality.
%
% Performs several classification tasks, such as determining directional
% tuning during eye and reach movements, determining behavioural state (ie.
% reaching, eye movement, planning). For BT M.Eng Thesis.
% 
% Use MATLAB's Naive Bayes Classifier implementation to classify eye
% movements into a category defined by type (saccade or pursuit)
% and direction (right, up, left, down).


fprintf('Beginning classification2_3.m!\n');
clear all;

%%%%%%%%%% TO SET %%%%%%%%%%%
% Choose which channels to look at
mipChans = 1:36; pmdChans = 37:48; allChans = 1:48;
channels = mipChans;
chanType = 'mip';%'pmd','both'
numChannels = length(channels);
pVals = [.6 .04 .004 .001 .0007];
numPs = length(pVals);

% Load an LFP data structure into workspace for trial info
%%%%%%%%%% TO SET %%%%%%%%%%%
type = 'sacs'; %purs,both
folder = 'M20100302_645';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Find number of, and particular, appropriate trials based on 'type'
% NOTE: THIS PART MUST SELECT SAME TRIALS AS anova_code4.m !!!
validTrials = 1:length(info);
if strcmp(type,'sacs')
    validTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'purs')
    validTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'both')
    validTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))) & info(validTrials,1)>128 & info(validTrials,1)<133);
end
numTrials = length(validTrials);

% Create labels vector that we can tweak later when we need to leave out
% certain trials
labels = zeros(numTrials,1);
for i = 1:numTrials
    if info(validTrials(i),1)==129 && info(validTrials(i),4)==65
        labels(i) = 1;
    elseif info(validTrials(i),1)==130 && info(validTrials(i),4)==65
        labels(i) = 2;
    elseif info(validTrials(i),1)==131 && info(validTrials(i),4)==65
        labels(i) = 3;
    elseif info(validTrials(i),1)==132 && info(validTrials(i),4)==65
        labels(i) = 4;
    elseif info(validTrials(i),1)==129 && info(validTrials(i),4)==66
        labels(i) = 5;
    elseif info(validTrials(i),1)==130 && info(validTrials(i),4)==66
        labels(i) = 6;
    elseif info(validTrials(i),1)==131 && info(validTrials(i),4)==66
        labels(i) = 7;
    elseif info(validTrials(i),1)==132 && info(validTrials(i),4)==66
        labels(i) = 8;
    end
end

% Load relevant channel spectrograms into workspace
% NOTE: spects matrix is cropped concatenation of spectrograms from all
% channels. The time dimension (3) begins at index 201 of the full
% spectrogram matrices, so to acces idx 281, use idx 81 in spect. (This
% note is only true for the original (mip and) pmdLookSpects.mat files).
if strcmp(chanType,'mip')
    %handle = sprintf('../data/%s/mipLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    load(handle);
elseif strcmp(chanType,'pmd')
    %handle = sprintf('../data/%s/pmdLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    load(handle);
else % Load relevant channels (all for now)
    %handle = sprintf('../data/%s/mipLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    mip=load(handle);
    %handle = sprintf('../data/%s/pmdLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    clear mip pmd;
end

% Load the tuning data obtained from anova_code4.m, which performed ANOVA
% for sacs and purs of each channel numTrials times while leavining out the 
% i-th trial each time.
handle = sprintf('../analysis/workspaces/%s/LOO_ANOVA_NN_160216.mat',folder);
tunData = load(handle);
% Make some helpful cell arrays for later
sigChans = cell(size(tunData.sacSigChans,2),size(tunData.sacSigChans,3));

% Choose which time steps we want to look. These are indices along the
% time dimension of the spects (later: spect, freqBands) matrix.
timeIdxs = [67 69 71 73 75]';
numTimeSteps = length(timeIdxs);
% Create the corresponding vector for indices of the tunData.times vector
tunTIdxs = zeros(numTimeSteps,1);
% How many freq bands will we look at
numBands = length(tunData.fBounds);

% Declare some matrices that will store output info to be saved
yTest = zeros(numTrials,1);
yOut = zeros(numTrials,numPs);

% For cases when using both look types, need to increment tunData.sac* and
% tunData.pur* separately
sTrial = 1; pTrial = 1;

% Enter leave-one-out loop
for l = 1:numTrials
    % Enter p-values loop
    for p = 1:numPs
        % Determine which channels are tuned for this pVal and (left out) trial
        for i = 1:numBands
            for j = 1:length(tunData.times)
                if strcmp(type,'sacs')
                    sigChans{i,j}=find(tunData.sacTuningPVals(channels,l,i,j)<pVals(p));
                elseif strcmp(type,'purs')
                    sigChans{i,j}=find(tunData.purTuningPVals(channels,l,i,j)<pVals(p));
                elseif strcmp(type,'both')
                    sigChans{i,j}=intersect(find(tunData.sacTuningPVals(channels,sTrial,i,j)<pVals(p)),find(tunData.purTuningPVals(channels,pTrial,i,j)<pVals(p)));
                    if p==numPs && i==numBands && j==length(tunData.times) && info(validTrials(l),4) == 65
                        if sTrial < size(tunData.sacTuningPVals,2)
                            sTrial = sTrial+1;
                        end
                    elseif p==numPs && i==numBands && j==length(tunData.times) && info(validTrials(l),4) == 66
                        if pTrial < size(tunData.purTuningPVals,2)
                            pTrial = pTrial+1;
                        end
                    end

                else
                    fprintf('Not a valid movement type!');
                end
            end
        end
        % Make sure there is something relevant in sigChans for pVal(p)
        if p > 1 && isempty(sigChans(:,tunTIdxs))
            fprintf('sigChans is empty for trial %d, pVal %d. Skipping...\n',l,p);
            continue;
        end
        
        % Leave out the l-th trial for testing
        trainTrials = [1:l-1 l+1:numTrials]';
        testTrials = l;
        
        % Begin building training matrix X, testing matrix Xtest
        
        % Declare X matrices to store samples to be classified
        % The last dimension is somewhat undetermined as we don't know how many
        % channels will be tuned in which frequency bands. Should scale okay if we
        % underestimate
        X = zeros(numTrials-1,250);
        Xtest = zeros(1,250);
        % Simply keep track of indices
        trainInd = 1; testInd = 1;
        % Loop through channels
        for chan = 1:numChannels
            % Get relevant part of spects
            spect = s(spects(chan,:,:,validTrials));
            % Separate spects into bands between 0-150 Hz, avoiding noise.
            [freqBands, fBounds] = create_freqBands(spect, f, tunData.fBounds);
            % On the first time through, find the indices in tunData.times
            % that match the times you want to use from freqBands (t). The
            % latter were chosen way up to in timeIdxs.
            if l==1 && p == 1
                for ts = 1:numTimeSteps
                    tunTIdxs(ts) = find(t(timeIdxs(ts)) == tunData.times);
                end
            end
            % If channels(chan) is tuned at time t and band f, use for training
            for ts = 1:numTimeSteps
                for fb = 1:numBands
                    % If current channel is among those considered "tuned", then
                    % use trials as training samples
                    if max(sigChans{fb,tunTIdxs(ts)}==chan) > 0
                        X(:,trainInd) = s(freqBands(fb,timeIdxs(ts),trainTrials));
                        trainInd=trainInd+1;
                        Xtest(:,testInd) = s(freqBands(fb,timeIdxs(ts),testTrials));
                        testInd=testInd+1;
                    end
                end
            end
            %fprintf('Done channel %d of %d\n',chan,numChannels);
        end
        % Adapt "answers" vectors from labels vector
        y = [labels(1:l-1); labels(l+1:end)];
        yTest(l) = labels(l);
        
        % Classify
        fprintf('Training Naive Bayes Classifier and Classifying!\n');
        
        % Train Naive Bayes model with X samples matrix
        NBModel =  class_train(X(:,1:trainInd-1),y);
        % Classify samples given in Xtest according to trained model and
        % keep track according to pVal(p) and left out trial l
        yOut(l,p) =  class_test(NBModel,Xtest(1:testInd-1));
        
        fprintf('Completed p-value %d of %d!\n',p,numPs);
    end
    fprintf('Completed classification for left out trial %d of %d!\n',l,numTrials);
end
answerCheck = zeros(numPs,numTrials,3);
numCorr = zeros(numPs,1); percentCorr = zeros(numPs,1);
answerCheck(:,:,1) = repmat(yTest,1,numPs)';
answerCheck(:,:,2) = yOut';
for i = 1:numPs
    answerCheck(i,:,3) = s(answerCheck(i,:,1))==s(answerCheck(i,:,2));
    numCorr(i) = s(sum(answerCheck(i,:,3)));
    percentCorr(i) = numCorr(i)/numTrials;
end

percentChance = 1/NBModel.NClasses;

%%%%%%%%%%%% TO SET %%%%%%%%%%%
descript = 'Time: 245, 300ms windows, all mip chans, all trials, sacs';
handle = sprintf('../analysis/workspaces/%s/mip_dirSacClass_LOO_170216_3.mat',folder);
save(handle,'descript','answerCheck','numCorr','percentCorr','percentChance','pVals','channels','type','validTrials','fBounds','t','timeIdxs');

% % Some analysis
figure;
plot(pVals,percentCorr);
xlabel('p-value');ylabel('Correct Prediction %');
title('Sac Direction Prediction Rate vs Confidence Level (MIP-LOO)');
xlim([0 .99]);
ax = gca; set(ax,'XScale','log');
ax.XTick = fliplr(pVals);
