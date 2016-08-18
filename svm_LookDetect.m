% svm_LookDetect.m
% Brendan Thorn - 02/16
%
% svm_LookDetect.m uses a support vector machine classifier to determine
% if, at a given time step, an eye movement (saccade or pursuit) will occur
% imminently.
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
handle = sprintf('../data/%s/mipLookSpects_100216.mat',folder);
%load(handle);
% mip=load(handle);
% handle = sprintf('../data/%s/pmdLookSpects',folder);
% pmd=load(handle);
% spects = [mip.spects; pmd.spects];
% t = mip.t;f=mip.f;
% clear mip pmd;

% Find number of, and particular, trials based on 'type' variables
validTrials = 1:length(info);

validSacTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133 & validTrials'~=302);
validPurTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
validBothTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))) & info(validTrials,1)>128 & info(validTrials,1)<133);

% Set bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15; 
bounds(3) = 25; bounds(4) = 35; 
bounds(5) = 45; bounds(6) = 55; 
bounds(7) = 65; bounds(8) = 75; 
bounds(9) = 85; bounds(10) = 95; 
bounds(11) = 105;bounds(12) = 125;
bounds(13) = 150;
numFreqBands = length(bounds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 020216 BT
% Try using SVM to classify eye movement vs not.
fprintf('Using SVM to detect eye movement times.\n');

%Get some details in order
numChans = size(spects,1);
type = 'sacs'; %%%%%%%%%%%%%%%%%%%%%%% <--- SET THIS 
if strcmp(type,'sacs')
    valTrials = validSacTrials;
elseif strcmp(type,'purs')
    valTrials = validPurTrials;
else
    valTrials = validBothTrials;
end

% TO SET
% Choose which trials to use (idxs within valTrials)
trials = 1:length(valTrials); numTrials = length(valTrials);
% Choose times to use for training and for testing
trainTimes = [10:2:30 32:39 42:43 48:64 66:2:84]';
testTimes = (9:2:99)';
%trainTimes = [28:3:82 86:91 98:101 102:3:180]';
%testTimes = (28:3:187)';
numTrains = length(trainTimes); numTests = length(testTimes);
% Note which training times should be classified as "MOVE"
moves = [42:43]';
% Note which testing times should be accepted as "MOVE"
testMoves = [41 43 45]';
% Choose some other parameters
polyOrds = [6]; numPolys = length(polyOrds);
boxCons = [.0001];% 10 100 1000 10000];
numBoxCons = length(boxCons);
% Limit the number of trials you want to work with
minTestTrial = 1; maxTestTrial = 100;
numTrials = maxTestTrial-minTestTrial+1;
% Save some results later
numTestTrials = 1; numTrainTrials = numTrials-numTestTrials; 
results = zeros(numPolys,numBoxCons,numTrials,numTests*numTestTrials,4);
numChans=24;
%SVMModels = cell(numPolys,numBoxCons,numTrials,numChans);

% Loop through trials to leave out for testing
for testTrial = minTestTrial:maxTestTrial
    % Set trainTrials and testTrials
    trainTrials = [valTrials(1:testTrial-1) valTrials(testTrial+1:numTrials)]';
    testTrials = valTrials(testTrial);
    numTrainTrials = length(trainTrials); numTestTrials = length(testTrials);
    % Zero svm matrices
    svmX = zeros(numTrains*numTrainTrials,numChans*numFreqBands);
    svmY = zeros(numTrains*numTrainTrials,1);
    testX = zeros(numTests*numTestTrials,numChans*numFreqBands);
    testY = zeros(numTests*numTestTrials,1);
    trialY = zeros(numTrains,1); testTrialY = zeros(numTests,1);
    trainRowInd = 1; testRowInd = 1; colInd = 1;
    % Loop through channels
    for chan = 1:numChans
        % reset rowInds
        trainRowInd = 1; testRowInd = 1;
        % Separate spectrogram into frequency bands
        spect = s(spects(chan,:,:,:));
        [freqBands, fBounds] = create_freqBands(spect, f, bounds);
        % Normalize each time step to same constant power
        %normBands = normByTime(freqBands);
        % Get average log powers at trainTimes for training trials
        svmX(:,colInd:colInd+numFreqBands-1) = reshape(permute(freqBands(:,trainTimes,trainTrials),[2 3 1]),numTrains*numTrainTrials,numFreqBands);
        %svmX(:,colInd+numFreqBands:colInd+2*numFreqBands-1) = reshape(permute(normBands(:,trainTimes,trainTrials),[2 3 1]),numTrains*numTrainTrials,numFreqBands);
        % Make vector of classes -- 1="MOVE", 0="NO MOVE"
        if chan == 1
            trialY(trainTimes>=min(moves) & trainTimes<=max(moves)) = 1;
            trialY(trainTimes<min(moves) | trainTimes>max(moves)) = 0;
            svmY = repmat(trialY,numTrainTrials,1);
        end
        % Get average log powers at testTimes for testing trials
        testX(:,colInd:colInd+numFreqBands-1) = reshape(permute(freqBands(:,testTimes,testTrials),[2 3 1]),numTests*numTestTrials,numFreqBands);
        %testX(:,colInd+numFreqBands:colInd+2*numFreqBands-1) = reshape(permute(normBands(:,testTimes,testTrials),[2 3 1]),numTests*numTestTrials,numFreqBands);
        if chan == 1
            testTrialY(testTimes>=min(testMoves) & testTimes<=max(testMoves)) = 1;
            testTrialY(testTimes<min(testMoves) | testTimes>max(testMoves)) = 0;
            testY = repmat(testTrialY,numTestTrials,1);
        end
        colInd = colInd + numFreqBands;
        %colInd = colInd + 2*numFreqBands;
    end
    
    fprintf('Training SVM for trial %d.\n',testTrial);
    for polyCount = 1:numPolys
        for boxCount = 1:numBoxCons
            % Train classifier
            SVMModel = fitcsvm(svmX,svmY,'Standardize',true,'Prior',struct('ClassNames',[0,1],'ClassProbs',[.99,.01]),'BoxConstraint',boxCons(boxCount),'KernelFunction','polynomial','KernelScale','auto','PolynomialOrder',polyOrds(polyCount),'OutlierFraction',0);
            % Fit posterior function for prediction
            SVMModel = fitPosterior(SVMModel,svmX,svmY);
            % Predict on testing set
            [testYout,score] = predict(SVMModel,testX);
            % Save results
            results(polyCount,boxCount,testTrial,:,1) = testY; results(polyCount,boxCount,testTrial,:,2)=testYout;results(polyCount,boxCount,testTrial,:,3)=(results(polyCount,boxCount,testTrial,:,1)==results(polyCount,boxCount,testTrial,:,2));results(polyCount,boxCount,testTrial,:,4)=score(:,2);
            % Clear some space before saving some info
            %SVMModel = compact(SVMModel);
            %SVMModels{polyCount,boxCount,testTrial} = SVMModel;

            fprintf('Done cross-val trial %d of %d for boxConstraint %d of %d for poly order %d of %d.\n',testTrial,maxTestTrial,boxCount,numBoxCons,polyCount,numPolys);
        end
    end
end

corrNoMoves = zeros(size(results,1),size(results,2));
corrMoves = zeros(size(results,1),size(results,2));
for i = 1:size(results,1)
    for j = 1:size(results,2)
        corrNoMoves(i,j) = (sum(sum(s(results(i,j,minTestTrial:maxTestTrial,1:15,3))))+sum(sum(s(results(i,j,minTestTrial:maxTestTrial,21:end,3)))))/(numTrials*15+(maxTestTrial-minTestTrial+1)*26);
        corrMoves(i,j) = sum(max(results(i,j,minTestTrial:maxTestTrial,17:19,3),[],4))/numTrials;
    end
end
% corrNoMoves = zeros(size(results,1),size(results,2));
% corrMoves = zeros(size(results,1),size(results,2));
% for i = 1:size(results,1)
%     for j = 1:size(results,2)
%         corrNoMoves(i,j) = (sum(sum(s(results(i,j,minTestTrial:maxTestTrial,1:19,3))))+sum(sum(s(results(i,j,minTestTrial:maxTestTrial,24:end,3)))))/((maxTestTrial-minTestTrial+1)*19+(maxTestTrial-minTestTrial+1)*31);
%         corrMoves(i,j) = sum(max(results(i,j,minTestTrial:maxTestTrial,21:22,3),[],4))/(maxTestTrial-minTestTrial+1);
%     end
% end

fprintf('Completed SVM detection.\n');

res1=s(results(1,1,minTestTrial:maxTestTrial,:,1));
res2=s(results(1,1,minTestTrial:maxTestTrial,:,2));
res3=s(results(1,1,minTestTrial:maxTestTrial,:,3));
res4=s(results(1,1,minTestTrial:maxTestTrial,:,4));
%corrNoMove = (sum(sum(res3(:,1:19)))+sum(sum(res3(:,24:end))))/(size(res3,1)*19+size(res3,1)*31);
%corrMove = sum(max(res3(:,21:22),[],2))/size(res3,1);

% Make some plots

% Plot MOVE probability heat map for trials vs time
figure;imagesc(t(testTimes),1:numTestTrials,res4);colorbar();ax=gca;ax.YDir='normal';
xlabel('Time Relative to Saccade');ylabel('Test Trial');
title('Saccade Probability -- SVM Saccade Detection');line([0 0],[minTestTrial maxTestTrial],'Color','r','LineWidth',2);

figure;imagesc(t(testTimes),1:numTestTrials,res2);colorbar();ax=gca;ax.YDir='normal';
xlabel('Time Relative to Saccade');ylabel('Test Trial');
title('Saccade Detects -- SVM Saccade Detection');line([0 0],[minTestTrial maxTestTrial],'Color','r','LineWidth',2);

