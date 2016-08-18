%%Brendan Thorn - 02/16
% classification2_4.m
%
% 20/02/16 - This is a refactored version of classification2_3.m. It uses
% the same LOO ANOVA data, but instead of training with every trial, as of
% the time I am writing this it uses 50 randomly selected trials. The goal
% is no longer to produce a correct percentage vs p-value plot, but instead
% compare areas MIP and PMD. There will be some other minor changes as well. 
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


fprintf('Beginning classification2_4.m!\n');
clearvars -except spects t f;

%%%%%%%%%% TO SET %%%%%%%%%%%
% Choose which channels to look at. You can make custom combination of
% channels later, just make sure they are contained in the ones chosen now.
mipChans = 1:36; pmdChans = 37:48; allChans = 1:48;
valChannels = pmdChans;
chanType = 'pmd';%'pmd','both'
numValChannels = length(valChannels);
pVals = [.01];
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
% TO SET:
% Choose a (random) subset of validTrials for testing if you'd like
numValTestTrials = 100;
valTestIdxs = randperm(numTrials,numValTestTrials);
valTestTrials = validTrials(valTestIdxs);

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
    else
        labels(i) = 0;
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

% To SET:
% Choose which time steps we want to look. These are indices along the
% time dimension of the spects (later: spect, freqBands) matrix.
timeIdxs = [71 73 75]';
numTimeSteps = length(timeIdxs);
% Create the corresponding vector for indices of the tunData.times vector
tunTIdxs = zeros(numTimeSteps,1);
% How many freq bands will we look at
numBands = length(tunData.fBounds);

% For cases when using both look types, need to increment tunData.sac* and
% tunData.pur* separately
sTrial = 1; pTrial = 1;

% TO SET !!!!
% Choose a number of channels from chanType to use. The numChannels "most
% tuned" will be used later.
numChannels = 4;
% Either choose channels for whole run here, or more granularly below
%channels = channels(randperm(numValChannels,numChannels));
% Choose how many randomly selected trials you want to use for training
trainNums = 10:10:170; numTrains = length(trainNums);
% Choose how many iterations to run for each number of train trials for
% each test trial
numIts = 10;

% Declare some matrices that will store output info to be saved
yTest = zeros(numIts,numTrains,numValTestTrials);
yOut = zeros(numIts,numTrains,numValTestTrials);

for itNum = 1:numIts
    % Enter leave-one-out loop
    for l = 1:numValTestTrials
        % Choose testTrials for this iteration from valTestTrials
        testTrials = valTestTrials(l);
        testIdxs = valTestIdxs(l);
        numTestTrials = length(testTrials);
        % Enter p-values loop
        for p = 1:numPs
            % Determine which channels are tuned for this pVal and (left out) trial
            for i = 1:numBands
                for j = 1:length(tunData.times)
                    if strcmp(type,'sacs')
                        sigChans{i,j}=valChannels(tunData.sacTuningPVals(valChannels,testIdxs,i,j)<pVals(p));
                    elseif strcmp(type,'purs')
                        sigChans{i,j}=valChannels(tunData.purTuningPVals(valChannels,testIdxs,i,j)<pVals(p));
                    elseif strcmp(type,'both')
                        % Don't use this until you figure out how to handle
                        % test trials for both look types simultaneously
                        %sigChans{i,j}=intersect(valChannels(tunData.sacTuningPVals(valChannels,sTrial,i,j)<pVals(p)),valChannels(tunData.purTuningPVals(valChannels,pTrial,i,j)<pVals(p)));
                        %if p==numPs && i==numBands && j==length(tunData.times) && info(validTrials(l),4) == 65
                        %    if sTrial < size(tunData.sacTuningPVals,2)
                        %        sTrial = sTrial+1;
                        %    end
                        %elseif p==numPs && i==numBands && j==length(tunData.times) && info(validTrials(l),4) == 66
                        %    if pTrial < size(tunData.purTuningPVals,2)
                        %        pTrial = pTrial+1;
                        %    end
                        %end
                    else
                        fprintf('Not a valid movement type!');
                    end
                end
            end
            % On the first time through, find the indices in tunData.times
            % that match the times you want to use from freqBands (t). The
            % latter were chosen way up to in timeIdxs.
            if itNum == 1 && l==1 && p == 1
                for ts = 1:numTimeSteps
                    tunTIdxs(ts) = find(t(timeIdxs(ts)) == tunData.times);
                end
            end
            % Select the numChannels most tuned channels at the relevant
            % times according to the sigChans matrix. "Most tuned" is taken
            % to mean, "appears most often in sigChans during time steps we
            % are using".
            tunedChanCount = zeros(numValChannels,1);
            for chan = 1:numValChannels
                for ts = 1:numTimeSteps
                    for fb = 1:numBands
                        % If current channel is "tuned", add to tally
                        if max(sigChans{fb,tunTIdxs(ts)}==valChannels(chan)) > 0
                            tunedChanCount(chan,1) = tunedChanCount(chan,1)+1;
                        end
                    end
                end
            end
            [S,I] = sort(tunedChanCount,'descend');
            channels = I(1:numChannels);
            for nTCount = 1:numTrains
                % Choose how main trials to use for training
                numTrainTrials = trainNums(nTCount);
                % Choose training trials. Make sure you don't use the
                % trial you are going to test as one of the training trials.
                valTrainIdxs = find(validTrials~=testTrials);
                valTrainTrials = validTrials(valTrainIdxs);
                numValTrainTrials = length(valTrainTrials);
                trainIdxs = valTrainIdxs(randperm(numValTrainTrials,numTrainTrials));
                trainTrials = validTrials(trainIdxs);
                
                % If you are using a small numTrainTrials, make sure you
                % have at least 2 examples of each look type. If you have a
                % larger value then don't waste the time
                if numTrainTrials < 46
                    % Use a product of the direction strobe and move-type
                    % strobe as signatures for possible look classes. Find
                    % out which combinations exist, then make sure we have
                    % enough of each type.
                    lookClasses = unique(info(validTrials,1).*info(validTrials,4));
                    flag = 1;
                    while flag == 1
                        for i = 1:length(lookClasses)
                            if sum(info(trainTrials,1).*info(trainTrials,4)==lookClasses(i)) < 2
                                % If there aren't enough, resample trainTrials
                                trainIdxs = valTrainIdxs(randperm(numValTrainTrials,numTrainTrials));
                                trainTrials = validTrials(trainIdxs);
                                % Check all lookClasses again
                                flag = 1;
                                break;
                            end
                            flag = 0;
                        end
                    end
                end
                
                % Begin building training matrix X, testing matrix Xtest
                
                % Declare X matrices to store samples to be classified
                % The last dimension is somewhat undetermined as we don't know how many
                % channels will be tuned in which frequency bands. Should scale okay if we
                % underestimate
                X = zeros(numTrainTrials,250);
                Xtest = zeros(1,250);
                % Simply keep track of indices
                trainInd = 1; testInd = 1;
                % Loop through channels
                for chan = 1:numChannels
                    % Get relevant part of spects
                    spect = s(spects(channels(chan),:,:,validTrials));
                    % Separate spects into bands between 0-150 Hz, avoiding noise.
                    [freqBands, fBounds] = create_freqBands(spect, f, tunData.fBounds);
                    % If channels(chan) is tuned at time t and band f, use for training
                    for ts = 1:numTimeSteps
                        for fb = 1:numBands
                            % If current channel is among those considered "tuned", then
                            % use trials as training samples
                            if max(sigChans{fb,tunTIdxs(ts)}==valChannels(channels(chan))) > 0
                                X(:,trainInd) = s(freqBands(fb,timeIdxs(ts),trainIdxs));
                                trainInd=trainInd+1;
                                Xtest(:,testInd) = s(freqBands(fb,timeIdxs(ts),testIdxs));
                                testInd=testInd+1;
                            end
                        end
                    end
                    %fprintf('Done channel %d of %d\n',chan,numChannels);
                end
                % Adapt "answers" vectors from labels vector
                y = labels(trainIdxs);
                yTest(itNum,nTCount,l) = labels(testIdxs);
                
                % Classify
                %fprintf('Training Naive Bayes Classifier and Classifying!\n');
                
                % Train Naive Bayes model with X samples matrix
                NBModel =  class_train(X(:,1:trainInd-1),y);
                % Classify samples given in Xtest according to trained model and
                % keep track according to pVal(p) and left out trial l
                yOut(itNum,nTCount,l) =  class_test(NBModel,Xtest(1:testInd-1));
                
                %fprintf('Completed p-value %d of %d!\n',p,numPs);
                fprintf('Completed classification for tSize %d of %d for test trial %d of %d for iteration %d of %d!\n',nTCount,numTrains,l,numValTestTrials,itNum,numIts);
            end
        end
    end
end

answerCheck = zeros(numIts,numTrains,numValTestTrials,3);
numCorr = zeros(numIts,numTrains); percentCorr = zeros(numIts,numTrains);

for i = 1:numIts
    for j = 1:numTrains
        answerCheck(i,j,:,1) = yTest(i,j,:);
        answerCheck(i,j,:,2) = yOut(i,j,:);
        answerCheck(i,j,:,3) = s(answerCheck(i,j,:,1))==s(answerCheck(i,j,:,2));
        numCorr(i,j) = s(sum(answerCheck(i,j,:,3)));
        percentCorr(i,j) = numCorr(i,j)/numValTestTrials;
    end
end
meanPercentCorr = mean(percentCorr,1);

percentChance = 1/NBModel.NClasses;

%%%%%%%%%%%% TO SET %%%%%%%%%%%
clear allChans ax chan Events freqBands info LFP spect spects tunData TT;
descript = 'Time: 1000, 300ms windows, 4 pmd chans, 10:10:170 train size, 100 trials, tIdxs 71,73,75, 10 its p01';
handle = sprintf('../analysis/workspaces/%s/w_pmd%d_dirSac_010316.mat',folder,numChannels);
save(handle);

% Some analysis
figure;
scatter(trainNums,meanPercentCorr);
xlabel('Number of Training Trials');ylabel('Correct Prediction %');
title('Sac Direction Prediction Rate vs Number of Training Trials (4 MIP-LOO)');
xlim([0 180]);ylim([.2 .5]);
