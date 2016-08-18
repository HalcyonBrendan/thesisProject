%%Brendan Thorn - 02/16
% classification2_4.m
%
% For BT M.Eng Thesis.
% 
% Use MATLAB's Naive Bayes Classifier implementation to classify eye
% movements into a category defined by type (saccade or pursuit)
% and direction (right, up, left, down).

fprintf('Beginning classification2_5.m!\n');
clearvars -except spects t f;
% Get a random number for identification
randNum = randi(100);

%%%%%%%%%% TO SET %%%%%%%%%%%
folder = 'M20100301_246';
dirType = 'trans';
chanType = 'mip';%'pmd','both'
type = 'sacs'; %purs,both
pVals = .05;
numPs = length(pVals);
% Choose a number of channels from chanType to use. The numChannels "most
% tuned" will be used later.
numChannels = 12;

% Load spectrograms
mipChans = 1:32; pmdChans = 33:48; allChans = 1:48;
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
    t = mip.t; f = mip.f;
    clear mip pmd;
    valChannels = allChans;
end
numValChannels = length(valChannels);

% Load an LFP data structure into workspace for trial info
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% Load the tuning data obtained from anova_code4.m, which performed ANOVA
% for sacs and purs of each channel numTrials times while leavining out the 
% i-th trial each time.
handle = sprintf('../analysis/workspaces/%s/%s_dirANOVA_240416.mat',folder,dirType);
tunData = load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Find number of, and particular, appropriate trials based on 'type'
% NOTE: THIS PART MUST SELECT SAME TRIALS AS anova_code4.m !!!
% Do eye movement direction at the same time
validTrials = 1:length(info);
eyeDir = zeros(size(info,1),1);
if strcmp(type,'sacs')
    validTrials = tunData.sacs;
    valCisTrials = tunData.cisSacs; valTransTrials = tunData.transSacs;
    eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==11 & info(:,4)==65) = 1;
    eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==-9 & info(:,4)==65) = 2;
    eyeDir(Events.EyeHandCoord(:,1)==-10 & Events.EyeHandCoord(:,2)==1 & info(:,4)==65) = 3;
    eyeDir(Events.EyeHandCoord(:,1)==10 & Events.EyeHandCoord(:,2)==1 & info(:,4)==65) = 4;
elseif strcmp(type,'purs')
    validTrials = tunData.purs;
    valCisTrials = tunData.cisPurs; valTransTrials = tunData.transPurs;
    eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==11 & info(:,4)==66) = 1;
    eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==-9 & info(:,4)==66) = 2;
    eyeDir(Events.EyeHandCoord(:,1)==-10 & Events.EyeHandCoord(:,2)==1 & info(:,4)==66) = 3;
    eyeDir(Events.EyeHandCoord(:,1)==10 & Events.EyeHandCoord(:,2)==1 & info(:,4)==66) = 4;
elseif strcmp(type,'both')
    validTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))));
    eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==11 & info(:,4)~=67) = 1;
    eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==-9 & info(:,4)~=67) = 2;
    eyeDir(Events.EyeHandCoord(:,1)==-10 & Events.EyeHandCoord(:,2)==1 & info(:,4)~=67) = 3;
    eyeDir(Events.EyeHandCoord(:,1)==10 & Events.EyeHandCoord(:,2)==1 & info(:,4)~=67) = 4;
end
numTrials = length(validTrials);
numValCisTrials = length(valCisTrials);
numValTransTrials = length(valTransTrials);
if strcmp(dirType,'cis')
    valDtTrials = valCisTrials;
    numValDtTrials = numValCisTrials;
elseif strcmp(dirType,'trans')
    valDtTrials = valTransTrials;
    numValDtTrials = numValTransTrials;
end

% Separate spects into bands between 0-150 Hz, avoiding noise.
[freqBands, fBounds] = create_freqBands(spects, f, tunData.fBounds);
clear spects;

% Choose a number of validTrials for testing if you'd like
numValTestTrials = min(numTrials,100);

% To SET:
% Choose which time steps we want to look. These are indices along the
% time dimension of the spects (later: spect, freqBands) matrix.
timeIdxs = tunData.timeIdxs;
numTimeSteps = size(timeIdxs,1);
numEpochs = size(timeIdxs,2);
epIdxs = find(ismember(tunData.timeIdxs,timeIdxs));
% Since first epIdx is from before eye move cue, use epIdxs(2)
%epIdxs(1) = epIdxs(2);
% How many freq bands will we look at
numBands = length(tunData.fBounds);

% For cases when using both look types, need to increment tunData.sac* and
% tunData.pur* separately
sTrial = 1; pTrial = 1;

trainNums = min(round(4*numValDtTrials/5),100);
numTrains = length(trainNums);
% Choose how many iterations to run for each number of train trials for each test trial
numIts = 500;

% Declare some matrices that will store output info to be saved
yTest = zeros(numEpochs,numIts,numTrains,numValTestTrials);
yOut = zeros(numEpochs,numIts,numTrains,numValTestTrials);
valTestIdxs = zeros(numValTestTrials,numIts);
valTestTrials = zeros(numValTestTrials,numIts);

for itNum = 1:numIts
    % Enter leave-one-out loop
    valTestIdxs(:,itNum) = randperm(numTrials,numValTestTrials);
    %valTestIdxs = 1:numValTestTrials;
    valTestTrials(:,itNum) = validTrials(valTestIdxs(:,itNum));
    for l = 1:numValTestTrials
        % Choose testTrials for this iteration from valTestTrials
        testTrials = valTestTrials(l,itNum);
        if ismember(testTrials,valDtTrials)
            testIdxs = find(ismember(valDtTrials,testTrials));
        else
            testIdxs = 1;
        end
        numTestTrials = length(testTrials);
        
        for nTCount = 1:numTrains
            % Choose how many trials to use for training
            numTrainTrials = trainNums(nTCount);
            % Choose training trials. Make sure you don't use the trial you are going to test as one of the training trials.
            valTrainIdxs = find(valDtTrials~=testTrials)';
            valTrainTrials = valDtTrials(valTrainIdxs)';
            numValTrainTrials = length(valTrainTrials);
            trainIdxs = valTrainIdxs(randperm(numValTrainTrials,numTrainTrials));
            trainTrials = valDtTrials(trainIdxs)';
            
            % If you are using a small numTrainTrials, make sure you
            % have at least 2 examples of each look type. If you have a
            % larger value then don't waste the time
            if numTrainTrials < 45
                % Use a product of the direction strobe and move-type strobe as signatures for possible look classes. Find
                % out which combinations exist, then make sure we have enough of each type.
                lookClasses = unique(eyeDir(validTrials).*info(validTrials,4));
                flag = 1;
                while flag == 1
                    for i = 1:length(lookClasses)
                        if sum(eyeDir(trainTrials).*info(trainTrials,4)==lookClasses(i)) < 2
                            % If there aren't enough, resample trainTrials
                            trainIdxs = valTrainIdxs(randperm(numValTrainTrials,numTrainTrials));
                            trainTrials = valDtTrials(trainIdxs);
                            % Check all lookClasses again
                            flag = 1;
                            break;
                        end
                        flag = 0;
                    end
                end
            end
            % Loop through epochs at which you want to decode
            for epNum = 1:numEpochs
                % Enter p-values loop
                for p = 1:numPs
                    % Determine which bands are tuned for this pVal and (left out) trial
                    if strcmp(type,'sacs')
                        sigBands = s(tunData.sacTuningPVals(valChannels,testIdxs,:,epIdxs(epNum))<pVals);
                    elseif strcmp(type,'purs')
                        sigBands = s(tunData.purTuningPVals(valChannels,testIdxs,:,epIdxs(epNum))<pVals);
                    end
                    
                    % Select the numChannels most tuned channels at the relevant times according to the sigChans matrix. "Most tuned" is taken
                    % to mean, "appears most often in sigChans during time steps we are using".
                    tunedChanCount = sum(sigBands,2);
                    [S,I] = sort(tunedChanCount,'descend');
                    channels = I(1:numChannels);

                    % Begin building training matrix X, testing matrix Xtest
                    
                    % Declare X matrices to store samples to be classified
                    % The last dimension is somewhat undetermined as we don't know how many
                    % channels will be tuned in which frequency bands. Should scale okay if we
                    % underestimate
                    X = zeros(numTrainTrials,100);
                    Xtest = zeros(1,100);
                    % Simply keep track of indices
                    trainInd = 1; testInd = 1;
                    % Loop through channels
                    for chan = 1:numChannels
                        numTuned =sum(sigBands(channels(chan),:));
                        X(:,trainInd:trainInd+numTuned-1) = s(freqBands(channels(chan),sigBands(channels(chan),:),timeIdxs(epNum),trainTrials))';
                        trainInd = trainInd+numTuned;
                        Xtest(:,testInd:testInd+numTuned-1) = s(freqBands(channels(chan),sigBands(channels(chan),:),timeIdxs(epNum),testTrials));
                        testInd = testInd+numTuned;
                    end
                    X=X(:,1:trainInd-1);
                    Xtest=Xtest(:,1:testInd-1);
                    
                    % Adapt "answers" vectors from labels vector
                    y = eyeDir(trainTrials);
                    yTest(epNum,itNum,nTCount,l) = eyeDir(testTrials);
                    
                    % Train Naive Bayes model with X samples matrix
                    NBModel =  class_train(X,y);
                    % Classify samples given in Xtest according to trained model and
                    % keep track according to pVal(p) and left out trial l
                    yOut(epNum,itNum,nTCount,l) =  class_test(NBModel,Xtest);
                    
                    %fprintf('Completed p-value %d of %d!\n',p,numPs);
                    %fprintf('Completed classification for epoch %d of %d for test trial %d of %d for iteration %d of %d!\n',epNum,numEpochs,l,numValTestTrials,itNum,numIts);
                end
            end
        end
    end
    fprintf('Completed iteration %d of %d\n',itNum,numIts);
end

answerCheck = zeros(numEpochs,numIts,numTrains,numValTestTrials,3);
numCorr = zeros(numEpochs,numIts,numTrains); percentCorr = zeros(numEpochs,numIts,numTrains);

for e = 1:numEpochs
    for i = 1:numIts
        for j = 1:numTrains
            answerCheck(e,i,j,:,1) = yTest(e,i,j,:);
            answerCheck(e,i,j,:,2) = yOut(e,i,j,:);
            answerCheck(e,i,j,:,3) = s(answerCheck(e,i,j,:,1))==s(answerCheck(e,i,j,:,2));
            numCorr(e,i,j) = s(sum(answerCheck(e,i,j,:,3)));
            percentCorr(e,i,j) = numCorr(e,i,j)/numValTestTrials;
        end
    end
end
meanPercentCorr = mean(percentCorr,2);

percentChance = 1/NBModel.NClasses;

%%%%%%%%%%%% TO SET %%%%%%%%%%%
handle = sprintf('../analysis/workspaces/%s/%s%s_dirClass_%dchs_240416_%s.mat',folder,type,chanType,numChannels,dirType);
save(handle,'dirType','answerCheck','channels','chanType','epIdxs','eyeDir','f','fBounds','folder','handle','I','trainNums','meanPercentCorr','mipChans','numCorr','percentCorr','pmdChans','pVals','S','sigBands','t','timeIdxs','trainIdxs','trainTrials','tunedChanCount','type','valChannels','validTrials','valTestIdxs','valTestTrials','valTrainIdxs','valTrainTrials','yOut','yTest','info','valCisTrials','valTransTrials');

% Some analysis
figure;
plot(t(timeIdxs),s(meanPercentCorr));
xlabel('Time');ylabel('Correct Prediction %');title(numChannels);
%title('Sac Direction Prediction Rate vs Number of Training Trials (4 MIP-LOO)');
%xlim([0 180]);ylim([.2 .5]);
