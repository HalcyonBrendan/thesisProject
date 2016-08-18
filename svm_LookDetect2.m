% svm_LookDetect2.m
% Brendan Thorn - 02/16
%
% svm_LookDetect2.m uses a support vector machine classifier to determine
% if, at a given time step, an eye movement (saccade or pursuit) will occur
% imminently. It is an adapted version of svm_LookDetect.m that is meant to
% be compatible with anovaDetect.m/manovaDetct.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Beginning svm_LookDetect.m at time:\n');
c=clock; fprintf('%.0f:%.0f\n',c(4),c(5));

fprintf('Loading data into workspace ...\n');
% Key variables to load before getting started
folder = 'M20100302_645';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
info = Events.dataSummary;
%%%%%%%%%% TO SET %%%%%%%%%%%
% Choose which channels to look at. You can make custom combination of
% channels later, just make sure they are contained in the ones chosen now.
mipChans = 1:32; pmdChans = 33:48; allChans = 1:48;
chanType = 'mip';%'pmd','both'
if strcmp(chanType,'mip')
    handle = sprintf('../data/%s/mipLookSpects_240_024.mat',folder);
    load(handle);
    valChannels = mipChans;
elseif strcmp(chanType,'pmd')
    handle = sprintf('../data/%s/pmdLookSpects_240_024.mat',folder);
    load(handle);
    valChannels = pmdChans;
else % Load relevant channels (all for now)
    handle = sprintf('../data/%s/mipLookSpects_240_024.mat',folder);
    mip=load(handle);
    handle = sprintf('../data/%s/pmdLookSpects_240_024.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    clear mip pmd;
    valChannels = allChans;
end
numValChans = length(valChannels);

% Load ANOVA information from anovaDetect.m
handle = sprintf('../analysis/workspaces/%s/LOO_detANOVA_NN_240_170316_6162.mat',folder);
detAnova = load(handle);
sigChans = cell(size(detAnova.sacSigChans,2),1);

% Find number of, and particular, trials based on 'type' variables
validTrials = 1:length(info);

% Set bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15; 
bounds(3) = 25; bounds(4) = 35; 
bounds(5) = 45; bounds(6) = 55; 
bounds(7) = 65; bounds(8) = 75; 
bounds(9) = 85; bounds(10) = 95; 
bounds(11) = 105;bounds(12) = 125;
bounds(13) = 150;
numFreqBands = length(bounds);

%Get some details in order
type = 'purs'; %%%%%%%%%%%%%%%%%%%%%%% <--- SET THIS 
if strcmp(type,'sacs')
    valTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133 & validTrials'~=302);
elseif strcmp(type,'purs')
    valTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
else
    valTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))) & info(validTrials,1)>128 & info(validTrials,1)<133 & validTrials'~=302);
end
numValTrials = length(valTrials);

% TO SET
% Choose a number of channels from chanType to use. The numChannels "most tuned" will be used later.
numChans = 12;
% Note which times (idxs) are classified as "MOVE" in ANOVA
moves = detAnova.moveIdxs;
detTIdxs = detAnova.tIdxs';
% Choose times (idxs) from ANOVA data to use for training and for testing
%trainTimes = [detTIdxs(1):moves(1)-3 moves(1):moves(end)+1 moves(end)+4:detTIdxs(end)]';
trainTimes = detTIdxs;
testTimes = 21:60;
numTrainTimes = length(trainTimes); numTestTimes = length(testTimes);
% Note which training/testing idxs should be accepted as "MOVE"
trainMoves = moves;
testMoves = 56:60;
% Determine which indices in trainTimes and testTimes that trainMoves and testMoves correspond to
trainMovesIdxs = find(ismember(trainTimes,trainMoves));
testMovesIdxs = find(ismember(testTimes,testMoves));
% Choose some other parameters
polyOrds = [7 8]; numPolys = length(polyOrds);
boxCons = [.003 .005 .007];numBoxCons = length(boxCons);
pVals = [.0001]; numPs = length(pVals);
classProbs = [.99 .01; .95 .05; .9 .1; .85 .15; .8 .2]; numCPs = size(classProbs,1);
outFracs = [0]; numOFs = length(outFracs);
% Limit the number of trials you want to work with
minTrainIdx = 1; maxTrainIdx = numValTrials;
trainIdxs = minTrainIdx:maxTrainIdx;
trainTrials = valTrials(trainIdxs); numTrainTrials = length(trainTrials);
minTestIdx = 1; maxTestIdx = length(valTrials);

% TO SET !!!
% Decide how many trials to test and train on at once
numCurrTestTrials = 1; numCurrTrainTrials = min(100,numTrainTrials);
% Choose which trials to test on
testIdxs = maxTestIdx-129:maxTestIdx-80; %%% <-- SET THIS
testTrials = valTrials(testIdxs); numTestTrials = length(testTrials);
% Store binary outcomes (MOVE or NOMOVE) in results matrix
results = zeros(numTestTrials,numPs,numPolys,numBoxCons,numCPs,numOFs,numTestTimes*numCurrTestTrials);
% Store outcome probabilities in scores matrix
scores = zeros(numTestTrials,numPs,numPolys,numBoxCons,numCPs,numOFs,numTestTimes*numCurrTestTrials);

% Loop through trials to leave out for testing
for testCount = 1:numTestTrials
    
    %Enter p-values loop
    for p = 1:numPs
        % Determine which channels are tuned for this pVal and (left out) trial
        for i = 1:numFreqBands
            if strcmp(type,'sacs')
                sigChans{i,1}=find(detAnova.sacDetectPVals(valChannels,testIdxs(testCount),i)<pVals(p));
            elseif strcmp(type,'purs')
                sigChans{i,1}=find(detAnova.purDetectPVals(valChannels,testIdxs(testCount),i)<pVals(p));
            elseif strcmp(type,'both')
                % TODO: implement the sTrial/pTrial incrementation properly if needed
                sigChans{i,1}=intersect(find(detAnova.sacDetectPVals(valChannels,sTrial,i)<pVals(p)),find(detAnova.purDetectPVals(valChannels,pTrial,i)<pVals(p)));
            end
        end
        % Select the numChannels most tuned channels at the relevant times according to the sigChans matrix. "Most tuned" is taken
        % to mean, "appears most often in sigChans during time steps we are using".
        tunedChanCount = zeros(numValChans,1);
        for chan = 1:numValChans
            for fb = 1:numFreqBands
                % If current channel is "tuned", add to tally
                if max(sigChans{fb,1}==valChannels(chan)) > 0
                    tunedChanCount(chan,1) = tunedChanCount(chan,1)+1;
                end
            end
        end
        [S,I] = sort(tunedChanCount,'descend');
        channels = I(1:numChans);
        % Zero svm matrices (it's okay if they are too big for now)
        trainX = zeros(numTrainTimes*numTrainTrials,numChans*numFreqBands);
        testX = zeros(numTestTimes*numTrainTrials,numChans*numFreqBands);
        trainTrialY = zeros(numTrainTimes,1); testTrialY = zeros(numTestTimes,1);
        rowInd = 1;colInd = 1;
        % Loop through channels
        for chan = 1:numChans
            % Reset rowInd
            rowInd = 1;
            % Determine for which frequency bands channels(chan) is tuned
            tunedBands = zeros(numFreqBands,1);
            for fb = 1:numFreqBands
                if sum(channels(chan)==sigChans{fb,1})>0
                    tunedBands(fb,1) = 1;
                end
            end
            numTunedBands = sum(tunedBands);
            % If no bands are tuned for this channel, go to next one
            if numTunedBands == 0
                continue;
            end
            tunedBands = find(tunedBands);

            % Separate spectrogram into frequency bands
            spect = s(spects(channels(chan),:,:,:));
            [freqBands, fBounds] = create_freqBands(spect, f, bounds);
            % Get average log powers at trainTimes for training trials
            trainX(:,colInd:colInd+numTunedBands-1) = reshape(permute(freqBands(tunedBands,trainTimes,valTrials),[2 3 1]),numTrainTimes*numTrainTrials,numTunedBands);
            testX(:,colInd:colInd+numTunedBands-1) = reshape(permute(freqBands(tunedBands,testTimes,valTrials),[2 3 1]),numTestTimes*numTrainTrials,numTunedBands);
            % Make vector of classes -- 1="MOVE", 0="NO MOVE"
            if chan == 1
                trainTrialY(trainTimes>=min(trainMoves) & trainTimes<=max(trainMoves)) = 1;
                trainY = repmat(trainTrialY,numValTrials,1);
                testTrialY(testTimes>=min(testMoves) & testTimes<=max(testMoves)) = 1;
                testY = repmat(testTrialY,numValTrials,1);
            end
            colInd = colInd + numTunedBands;
        end
        trainX = trainX(:,1:colInd-1);
        testX = testX(:,1:colInd-1);
        
        % Get the correct rows from the trainX/testX matrices for training
        % and testing given the current test trial
        notX = trainX(1+(testCount-1)*numTrainTimes:testCount*numTrainTimes,:);
        [X,idxs] = setdiff(trainX,notX,'rows');
        Y = trainY(idxs);
        tX = testX(1+(testCount-1)*numTestTimes:testCount*numTestTimes,:);
        tY = testY(1+(testCount-1)*numTestTimes:testCount*numTestTimes,1);
        
        %fprintf('Running SVM with parameter variations for trial %d.\n',testCount);
        for pCs = 1:numPolys
            for bCs = 1:numBoxCons
                for cPs = 1:numCPs
                    for oFs = 1:numOFs
                        % Train classifier
                        SVMModel = fitcsvm(X,Y,'Standardize',true,'Prior',struct('ClassNames',[0,1],'ClassProbs',classProbs(cPs,:)),'BoxConstraint',boxCons(bCs),'KernelFunction','polynomial','KernelScale','auto','PolynomialOrder',polyOrds(pCs),'OutlierFraction',outFracs(oFs),'IterationLimit',2e4,'CacheSize','maximal');
                        % Fit posterior function for prediction
                        SVMModel = fitPosterior(SVMModel,X,Y);
                        % Predict on testing set and save results
                        [result,score] = predict(SVMModel,tX);
                        results(testCount,p,pCs,bCs,cPs,oFs,:) = result;
                        scores(testCount,p,pCs,bCs,cPs,oFs,:) = s(score(:,2));

                        %fprintf('Num support vectors: %d of %d total.\n',size(SVMModel.SupportVectors,1),size(SVMModel.X,1));
                    end
                end
            end
            fprintf('Done cross-val trial %d of %d for poly order %d of %d for p-value %d of %d.\n',testCount,numTestTrials,pCs,numPolys,p,numPs);
        end
        fprintf('Done cross-val trial %d of %d for p-value %d of %d.\n',testCount,numTestTrials,p,numPs);
    end
end

% Now loop through all of the parameters again and get some useful
% statistics regarding the SVM performance for each parameter set
truePos = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs);trueNegs = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs);falsePos = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs);falseNegs = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs);MCC = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs);
movesCorrPerc = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs); noMovesCorrPerc = zeros(numPs,numPolys,numBoxCons,numCPs,numOFs);
for p = 1:numPs
    for pCs = 1:numPolys
        for bCs = 1:numBoxCons
            for cPs = 1:numCPs
                for oFs = 1:numOFs
                    % Get number of true positives (defined as a correct detection at AT LEAST one of the times given in testMoves)
                    truePos(p,pCs,bCs,cPs,oFs) = sum(max(s(results(:,p,pCs,bCs,cPs,oFs,testMovesIdxs)),[],2));
                    % Get number of true negatives (defined as a correct non-detection at any tested time not given in testMoves)
                    trueNegs(p,pCs,bCs,cPs,oFs) = sum(sum(s(results(:,p,pCs,bCs,cPs,oFs,1:testMovesIdxs(1)-1))==0))+sum(sum(s(results(:,p,pCs,bCs,cPs,oFs,testMovesIdxs(end)+1:end))==0));
                    % Get number of false negatives (defined as a missed detection at ALL of the times given in testMoves)
                    falseNegs(p,pCs,bCs,cPs,oFs) = numTestTrials - truePos(p,pCs,bCs,cPs,oFs);
                    % Get number of false positives (defined as an incorrect detection at any tested time not given in testMoves)
                    falsePos(p,pCs,bCs,cPs,oFs) = sum(sum(s(results(:,p,pCs,bCs,cPs,oFs,1:testMovesIdxs(1)-1))))+sum(sum(s(results(:,p,pCs,bCs,cPs,oFs,testMovesIdxs(end)+1:end))));
                    % Compute Matthews correlation coefficient for this parameterization of the SVM
                    MCC(p,pCs,bCs,cPs,oFs) = computeMatthews(truePos(p,pCs,bCs,cPs,oFs),trueNegs(p,pCs,bCs,cPs,oFs),falseNegs(p,pCs,bCs,cPs,oFs),falsePos(p,pCs,bCs,cPs,oFs));
                    
                    % Print out some results for fun
                    movesCorrPerc(p,pCs,bCs,cPs,oFs) = truePos(p,pCs,bCs,cPs,oFs)/(truePos(p,pCs,bCs,cPs,oFs)+falseNegs(p,pCs,bCs,cPs,oFs));
                    noMovesCorrPerc(p,pCs,bCs,cPs,oFs) = trueNegs(p,pCs,bCs,cPs,oFs)/(trueNegs(p,pCs,bCs,cPs,oFs)+falsePos(p,pCs,bCs,cPs,oFs));
                    fprintf('Results (%d,%d,%d,%d,%d): Correct Detect Percent: %f; Correct Non-Detect Percent: %f; MCC: %f;\n',p,pCs,bCs,cPs,oFs,movesCorrPerc(p,pCs,bCs,cPs,oFs),noMovesCorrPerc(p,pCs,bCs,cPs,oFs),MCC(p,pCs,bCs,cPs,oFs));
                end
            end
        end
    end
end

fprintf('Completed SVM detection.\n');

% Put some of the results into a more palatable format
res=s(results(:,1,1,1,1,1,:));
sco=s(scores(:,1,1,1,1,1,:));


% Make some plots

% Plot MOVE probability heat map for trials vs time
figure;imagesc(t(testTimes),1:numTestTrials,res);colorbar();ax=gca;ax.YDir='normal';
xlabel('Time Relative to Saccade');ylabel('Test Trial');
title('Saccade Probability -- SVM Saccade Detection');line([-.12 -.12],[0 numTestTrials],'Color','r','LineWidth',2);line([0 0],[0 numTestTrials],'Color','g','LineWidth',2);line([-.17 -.17],[0 numTestTrials],'Color','c','LineWidth',2);line([-.07 -.07],[0 numTestTrials],'Color','c','LineWidth',2);

figure;imagesc(t(testTimes),1:numTestTrials,sco);colorbar();ax=gca;ax.YDir='normal';
xlabel('Time Relative to Saccade');ylabel('Test Trial');
title('Saccade Detects -- SVM Saccade Detection');line([-.12 -.12],[0 numTestTrials],'Color','r','LineWidth',2);line([0 0],[0 numTestTrials],'Color','g','LineWidth',2);line([-.17 -.17],[0 numTestTrials],'Color','c','LineWidth',2);line([-.07 -.07],[0 numTestTrials],'Color','c','LineWidth',2);

