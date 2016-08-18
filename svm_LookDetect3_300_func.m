function [ ] = svm_LookDetect3_300_func( folder,numChans,chType,lookType,pVal,polyOrd,boxCon,classProb,randNum )
% svm_LookDetect2.m
% Brendan Thorn - 02/16

% 22/03/16 - Converted svm_LookDetect2_300.m to this function so it can be
% called from svm_LookDetect2_300_paramOpt.m with the optimal parameters
%
% svm_LookDetect2.m uses a support vector machine classifier to determine
% if, at a given time step, an eye movement (saccade or pursuit) will occur
% imminently. It is an adapted version of svm_LookDetect.m that is meant to
% be compatible with anovaDetect.m/manovaDetct.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Beginning svm_LookDetect.m at time:\n');
c=clock; fprintf('%.0f:%.0f\n',c(4),c(5));
rng(1); % For reproducibility

fprintf('Loading data into workspace ...\n');
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
info = Events.dataSummary;
%%%%%%%%%% TO SET %%%%%%%%%%%
% Choose which channels to look at. You can make custom combination of
% channels later, just make sure they are contained in the ones chosen now.
mipChans = 1:32; pmdChans = 33:48; allChans = 1:48;
chanType = chType;%'mip','pmd','both'
if strcmp(chanType,'mip')
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    load(handle);
    valChannels = mipChans;
elseif strcmp(chanType,'pmd')
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    load(handle);
    valChannels = pmdChans;
else % Load relevant channels (all for now)
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    mip=load(handle);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    clear mip pmd;
    valChannels = allChans;
end
numValChans = length(valChannels);

% Load ANOVA information from anovaDetect.m
handle = sprintf('../analysis/workspaces/%s/detANOVA_170416_6973v7781.mat',folder);
detAnova = load(handle);
sigChans = cell(size(detAnova.sacSigChans,2),1);

% Set bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15; 
bounds(3) = 25; bounds(4) = 35; 
bounds(5) = 45; bounds(6) = 55; 
bounds(7) = 65; bounds(8) = 75; 
bounds(9) = 85; bounds(10) = 95; 
bounds(11) = 105;bounds(12) = 125;
bounds(13) = 150;
numFreqBands = length(bounds);
% Separate spectrogram into frequency bands
[freqBands, fBounds] = create_freqBands(spects, f, detAnova.fBounds);
clear spects;

%Get some details in order
type = lookType;
if strcmp(type,'sacs')
    %valTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    valTrials = detAnova.validSacs;
elseif strcmp(type,'purs')
    %valTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
    valTrials = detAnova.validPurs;
else
    %valTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))) & info(validTrials,1)>128 & info(validTrials,1)<133);
    valTrials = unique([detAnova.validSacs detAnova.validPurs]);
end
numValTrials = length(valTrials);

% TO SET
% Note which training/testing idxs should be accepted as "MOVE"
trainMoves = 76:79;
testMoves = 75:79;
% Choose times (idxs) from ANOVA data to use for training and for testing
trainTimes = [39:2:73 trainMoves]';
testTimes = 39:79;
numTrainTimes = length(trainTimes); numTestTimes = length(testTimes);

% Determine which indices in trainTimes and testTimes that trainMoves and testMoves correspond to
trainMovesIdxs = find(ismember(trainTimes,trainMoves));
testMovesIdxs = find(ismember(testTimes,testMoves));
% Choose some other parameters
polyOrds = polyOrd; numPolys = length(polyOrds);
boxCons = boxCon;numBoxCons = length(boxCons);
pVals = pVal; numPs = length(pVals);
classProbs = classProb;
% Limit the number of trials you want to work with
minTrainIdx = 1; maxTrainIdx = numValTrials;
valTrainIdxs = minTrainIdx:maxTrainIdx;

% TO SET !!!
% Decide how many trials to test at once
numCurrTestTrials = 1;
% Choose which trials to test on
minTestIdx = max(round(numValTrials/3),40)+1; maxTestIdx = numValTrials; %% <-- SET THIS
testIdxs = minTestIdx:maxTestIdx;
testTrials = valTrials(testIdxs); numTestTrials = length(testTrials);
% Store binary outcomes (MOVE or NOMOVE) in results matrix
results = zeros(numTestTrials,numTestTimes*numCurrTestTrials);
% Store outcome probabilities in scores matrix
scores = zeros(numTestTrials,numTestTimes*numCurrTestTrials);

% Loop through trials to leave out for testing
for testCount = 1:numTestTrials
    % Check if the current test trial is in the set of training trials
    trainIdxs = valTrainIdxs(valTrainIdxs~=testIdxs(testCount));
    trainTrials = valTrials(trainIdxs); numTrainTrials = length(trainTrials);
    %Enter p-values loop
    for p = 1:numPs
        % Determine which bands are tuned for this pVal and (left out) trial
        if strcmp(type,'sacs')
            sigBands = s(detAnova.sacDetectPVals(valChannels,testIdxs(testCount),:)<pVals(p));
        elseif strcmp(type,'purs')
            sigBands = s(detAnova.purDetectPVals(valChannels,testIdxs(testCount),:)<pVals(p));
        end
        % Select the numChannels most tuned channels at the relevant times according to the sigChans matrix. "Most tuned" is taken
        % to mean, "appears most often in sigChans during time steps we are using".
        tunedChanCount = sum(sigBands,2);
        [S,I] = sort(tunedChanCount,'descend');
        channels = I(1:numChans);
        
        % Zero svm matrices (it's okay if they are too big for now)
        trainX = zeros(numTrainTimes*numTrainTrials,numChans*numFreqBands);
        testX = zeros(numTestTimes*numCurrTestTrials,numChans*numFreqBands);
        trainTrialY = zeros(numTrainTimes,1); testTrialY = zeros(numTestTimes,1);
        colInd = 1;
        % Loop through channels
        for chan = 1:numChans
            numTuned =sum(sigBands(channels(chan),:));

            % Get average log powers at trainTimes for training trials
            trainX(:,colInd:colInd+numTuned-1) = reshape(permute(s(freqBands(channels(chan),sigBands(channels(chan),:),trainTimes,trainTrials)),[2 3 1]),numTrainTimes*numTrainTrials,numTuned);
            testX(:,colInd:colInd+numTuned-1) = reshape(permute(s(freqBands(channels(chan),sigBands(channels(chan),:),testTimes,testTrials(testCount))),[2 3 1]),numTestTimes*numCurrTestTrials,numTuned);
            % Make vector of classes -- 1="MOVE", 0="NO MOVE"
            if chan == 1
                trainTrialY(trainTimes>=min(trainMoves) & trainTimes<=max(trainMoves)) = 1;
                trainY = repmat(trainTrialY,numTrainTrials,1);
                testTrialY(testTimes>=min(testMoves) & testTimes<=max(testMoves)) = 1;
                testY = repmat(testTrialY,numCurrTestTrials,1);
            end
            colInd = colInd + numTuned;
        end
        trainX = trainX(:,1:colInd-1);
        testX = testX(:,1:colInd-1);
        
        %fprintf('Running SVM with parameter variations for trial %d.\n',testCount);
        for pCs = 1:numPolys
            for bCs = 1:numBoxCons
                % Train classifier
                %SVMModel = fitcsvm(trainX,trainY,'Standardize',true,'Prior',struct('ClassNames',[0,1],'ClassProbs',classProbs(cPs,:)),'BoxConstraint',boxCons(bCs),'KernelFunction','polynomial','KernelScale','auto','PolynomialOrder',polyOrds(pCs),'OutlierFraction',outFracs(oFs),'IterationLimit',2e4,'CacheSize','maximal');
                SVMModel = fitcsvm(trainX,trainY,'Standardize',true,'Prior',struct('ClassNames',[0,1],'ClassProbs',classProbs),'BoxConstraint',boxCons(bCs),'KernelFunction','rbf','KernelScale','auto','OutlierFraction',0,'IterationLimit',5e4);
                % Fit posterior function for prediction
                SVMModel = fitPosterior(SVMModel,trainX,trainY);
                % Predict on testing set and save results
                [result,score] = predict(SVMModel,testX);
                results(testCount,:) = result;
                scores(testCount,:) = s(score(:,2));
                
                %fprintf('Num support vectors: %d of %d total.\n',size(SVMModel.SupportVectors,1),size(SVMModel.X,1));
            end
            fprintf('Done cross-val trial %d of %d for poly order %d of %d for p-value %d of %d.\n',testCount,numTestTrials,pCs,numPolys,p,numPs);
        end
        %fprintf('Done cross-val trial %d of %d for p-value %d of %d.\n',testCount,numTestTrials,p,numPs);
    end
end

% Adjust SVM threshold to fix false positive rate at fpRate*100%
fpRate = .05;
sortedScores = sort(reshape(scores(:,1:testMovesIdxs(1)-1),[],1),'descend');
numFPs = round(fpRate*length(sortedScores));
threshValue = sortedScores(numFPs);
adjResults = scores>threshValue;

% Get number of true positives (defined as a correct detection at AT LEAST one of the times given in testMoves)
truePos = sum(max(adjResults(:,testMovesIdxs),[],2));
% Get number of true negatives (defined as a correct non-detection at any tested time not given in testMoves)
trueNegs = sum(sum(adjResults(:,1:testMovesIdxs(1)-1)==0));
% Get number of false negatives (defined as a missed detection at ALL of the times given in testMoves)
falseNegs = numTestTrials - truePos;
% Get number of false positives (defined as an incorrect detection at any tested time not given in testMoves)
falsePos = sum(sum(adjResults(:,1:testMovesIdxs(1)-1)));
% Compute Matthews correlation coefficient for this parameterization of the SVM
MCC = computeMatthews(truePos,trueNegs,falseNegs,falsePos);

% Print out some results for fun
movesCorrPerc = truePos/(truePos+falseNegs);
noMovesCorrPerc = trueNegs/(trueNegs+falsePos);
fprintf('Results: Correct Detect Percent: %f; Correct Non-Detect Percent: %f; MCC: %f;\n',movesCorrPerc,noMovesCorrPerc,MCC);

            
fprintf('Completed SVM detection.\n');

% Make some plots

% Plot MOVE probability heat map for trials vs time
% figure;imagesc(t(testTimes),1:numTestTrials,scores);colorbar();ax=gca;ax.YDir='normal';
% xlabel('Time Relative to Saccade [s]');ylabel('Test Trial');
% tit = sprintf('%s %s Probs SVM 170416',chType,lookType);
% title(tit);line([-.15 -.15],[0 numTestTrials],'Color','r','LineWidth',2);line([0 0],[0 numTestTrials],'Color','g','LineWidth',2);line([-.05 -.05],[0 numTestTrials],'Color','c','LineWidth',2);
% 
% figure;imagesc(t(testTimes),1:numTestTrials,adjResults);colorbar();ax=gca;ax.YDir='normal';
% xlabel('Time Relative to Pursuit [s]');ylabel('Test Trial');
% tit = sprintf('%s %s Detects SVM 170416',chType,lookType);
% title(tit);line([-.15 -.15],[0 numTestTrials],'Color','r','LineWidth',2);line([0 0],[0 numTestTrials],'Color','g','LineWidth',2);line([-.05 -.05],[0 numTestTrials],'Color','c','LineWidth',2);

clear allChans ax bCs chan colInd s detAnova Events fb freqBands i info LFP mipChans oFs p pCs pmdChans spect spects SVMModel testCount testY testX trainIdxs trainTrials trainY trainX TT validTrials;
handle = sprintf('../analysis/workspaces/%s/%s%s_optRun_%dchs_240416.mat',folder,type,chanType,numChans);
save(handle);

end