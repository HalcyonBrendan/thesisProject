%%Brendan Thorn - 01/16
% classification2_1.m
%
% This is a reorganized version of classification2.m that requires
% less reading from disk to reduce runtime. The ordering (and indexing) is
% much less intuitive, but is necessary for practicality.
%
%
% Performs several classification tasks, such as determining directional
% tuning during eye and reach movements, determining behavioural state (ie.
% reaching, eye movement, planning). For BT M.Eng Thesis.
% 

% Use MATLAB's Naive Bayes Classifier implementation to classify eye
% movements into a category defined by type (saccade or pursuit)
% and direction (right, up, left, down).

% Use last three time steps with no eye movement contribution as features
% for training classifier and then performing classification ie. eye
% movements start at t=0, and since spectrogram resolution is 25ms and
% windows 300ms, centred at current time, the last three time steps with no
% eye movement contributions are centred at t=-.150,-.175,-.200, in reverse
% order. In the case that a window overlaps with t=0, we shift backwards to
% avoid contamination.

fprintf('Beginning classification2_1.m!\n');
clear all;

% Choose which channels to look at
channels = 1:48;
numChannels = length(channels);

% TO SET
numRuns = 50;
pVals = [.99:-.05:.29 .29:-.04:.07 .06:-.01:.01 .009:-.001:.001 .0005];

numPs = length(pVals);

% Load an LFP data structure into workspace for trial info
%%%%%%%%%% TO SET %%%%%%%%%%%
type = 'both'; %purs,both
folder = 'M20100408_256';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Find number of, and particular, trials based on 'type' variables
validTrials = 1:length(info);
if strcmp(type,'sacs')
    validTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'purs')
    validTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'both')
    validTrials = validTrials(info(validTrials,4)~=67 & ~isnan(TT.T8(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
end
numTrials = length(validTrials);


% Run anova_code2 to determine the confidence level at which the channels
% are tuned at different times and frequencies
%anova_code2_func(pVals(1), 95, folder);
% Load tuning information produced by anova_code2.m
handle = sprintf('../analysis/workspaces/%s/tunedChans_%d_NN_jan0716.mat',folder,95);
tunChans = load(handle);

% Choose which time steps we want to look at for each band
% TODO: decide whether these are legitimate times to look at
timeIdxs = 281:2:289;
numTimeSteps = length(timeIdxs);
% How many freq bands will we look at
numBands = length(tunChans.fBounds);

% Before starting, a bit of code to put eye move times into one vector
alignTimes = zeros(numTrials,1);
for i = 1:numTrials
    if info(i,4)==65 & ~isnan(TT.T65sac)
        alignTimes(i)=TT.T65sac(i);
    elseif info(i,4)==66 & ~isnan(TT.T66on)
        alignTimes(i)=TT.T66on(i);
    else
        alignTimes(i)=TT.T8(i);
    end
end

% Prepare some constants/matrices that keep track of which trials
% were used for training/testing for each pVal and run, as we need
% to use the same ones each time we increment the channel
%%%%%%%%%%% TODO: Set up proper cross-validation %%%%%%%%%%%
% Choose approximate proportion of trials you want to train on (difficult
% to get exact number from here as we will later exclude fixation trials)
numTrainTrials = round(4/5*numTrials);
numTestTrials = numTrials-numTrainTrials;
trainingTrials = zeros(numPs,numRuns,numTrainTrials);
testingTrials = zeros(numPs,numRuns,numTestTrials);

% Declare X matrices to store samples to be classified
% The last dimension is somewhat undetermined as we don't know how many
% channels will be tuned in which frequency bands. Should scale okay if we
% underestimate
X = zeros(numPs,numRuns,numTrainTrials,250);
Xtest = zeros(numPs,numRuns,numTestTrials,250);

% Declare y answer vectors as well
y = zeros(numPs,numRuns,numTrainTrials);
yTest = zeros(numPs,numRuns,numTestTrials);
yOut = zeros(numPs,numRuns,numTestTrials);
answerCheck = zeros(numPs,numRuns,numTestTrials,3);
numCorr = zeros(numPs,numRuns);
percentCorr = zeros(numPs,numRuns);

% Keep track of index in X sample structures
trainInd = ones(numPs,numRuns); 
testInd = ones(numPs,numRuns);

for chan = 1:numChannels
    %%%%%%%%%%%%%%% TO SET %%%%%%%%%%%%%%%%%%
    % Load spectrogram for channels(chan)
    handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,channels(chan),47);
    % Get spectrogram/coherogram/cross-power for channel
    [spect,t,f] = getCohs(1,handle,channels(chan),0);
    % Center time vector
    t=t-7.5;
    fprintf('Spectrogram retrieved for channel %d.\n',channels(chan));
    
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
    % Or just pretend you normalized
    normSpect = log(permute(spect,[2 1 3]));
    clear spect;
    
    % Separate normalized spectrograms into bands between 0-200 Hz. Avoid noisy areas ie.
    % ~60 Hz and harmonics
    [freqBands, fBounds] = create_freqBands(normSpect, f, tunChans.fBounds);
    clear normSpect;
    
    % Loop through all desired pVals
    for pCounter = 1:numPs
        
        % Re-evaluate which channels are tuned for each movement type at
        % the level required be new pVal
        for i = 1:size(tunChans.sacTuningPVals,2)
            for j = 1:size(tunChans.sacTuningPVals,3)
                tunChans.sacSigChans{i,j}=find(tunChans.sacTuningPVals(channels,i,j)<pVals(pCounter));
                tunChans.purSigChans{i,j}=find(tunChans.purTuningPVals(channels,i,j)<pVals(pCounter));
            end
        end
        
        % Figure out which channels to use based on which are tuned for both
        % pursuits and saccades, or just choose one or the other
        if strcmp(type,'sacs')
            sigChans = tunChans.sacSigChans;
        elseif strcmp(type,'purs')
            sigChans = tunChans.purSigChans;
        elseif strcmp(type,'both')
            sigChans = cell(size(tunChans.purSigChans,1),size(tunChans.purSigChans,2));
            for i = 1:size(tunChans.purSigChans,1)
                for j = 1:size(tunChans.purSigChans,2)
                    sigChans{i,j} = intersect(tunChans.purSigChans{i,j},tunChans.sacSigChans{i,j});
                end
            end
        end
        
        
        % Repeat 
        for runCounter = 1:numRuns
            % For first channel set which trials will be used for what for
            % each run in each pVal
            if chan==1
                %%%%%%%%%%%%%%%%% TO SET %%%%%%%%%%%%%%%%%%%%
                % Randomly select trials that will be used for training
                randTrials = randperm(numTrials);
                trainingTrials(pCounter,runCounter,:) = validTrials(randTrials(1:numTrainTrials));
                % Trials for testing are those not used for training
                testingTrials(pCounter,runCounter,:) = validTrials(randTrials(numTrainTrials+1:end));

                % In case repeatability is desired, select first numTrainTrials
                %trainingTrials = validTrials(1:numTrainTrials);
                % In case repeatability is desired, select last numTestTrials
                %testingTrials = validTrials(numTrainTrials+1:numTrials);
            
            
                % Build "answers" vector for training samples
                % Each combination of up, down, left, right look directions, along with eye
                % movement type, saccade or pursuit, will be assigned a numerical value, as
                % follows:
                % Saccade-up     == 1
                % Saccade-down   == 2
                % Saccade-left   == 3
                % Saccade-right  == 4
                % Pursuit-up     == 5
                % Pursuit-down   == 6
                % Pursuit-left   == 7
                % Pursuit-right  == 8
                for i = 1:size(trainingTrials,3)
                    if info(trainingTrials(pCounter,runCounter,i),1)==129 && info(trainingTrials(pCounter,runCounter,i),4)==65
                        y(pCounter,runCounter,i) = 1;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==130 && info(trainingTrials(pCounter,runCounter,i),4)==65
                        y(pCounter,runCounter,i) = 2;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==131 && info(trainingTrials(pCounter,runCounter,i),4)==65
                        y(pCounter,runCounter,i) = 3;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==132 && info(trainingTrials(pCounter,runCounter,i),4)==65
                        y(pCounter,runCounter,i) = 4;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==129 && info(trainingTrials(pCounter,runCounter,i),4)==66
                        y(pCounter,runCounter,i) = 5;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==130 && info(trainingTrials(pCounter,runCounter,i),4)==66
                        y(pCounter,runCounter,i) = 6;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==131 && info(trainingTrials(pCounter,runCounter,i),4)==66
                        y(pCounter,runCounter,i) = 7;
                    elseif info(trainingTrials(pCounter,runCounter,i),1)==132 && info(trainingTrials(pCounter,runCounter,i),4)==66
                        y(pCounter,runCounter,i) = 8;
                    end
                end
                % Build "answers" vector for testing samples, as before
                for i = 1:size(testingTrials,3)
                    if info(testingTrials(pCounter,runCounter,i),1)==129 && info(testingTrials(pCounter,runCounter,i),4)==65
                        yTest(pCounter,runCounter,i) = 1;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==130 && info(testingTrials(pCounter,runCounter,i),4)==65
                        yTest(pCounter,runCounter,i) = 2;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==131 && info(testingTrials(pCounter,runCounter,i),4)==65
                        yTest(pCounter,runCounter,i) = 3;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==132 && info(testingTrials(pCounter,runCounter,i),4)==65
                        yTest(pCounter,runCounter,i) = 4;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==129 && info(testingTrials(pCounter,runCounter,i),4)==66
                        yTest(pCounter,runCounter,i) = 5;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==130 && info(testingTrials(pCounter,runCounter,i),4)==66
                        yTest(pCounter,runCounter,i) = 6;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==131 && info(testingTrials(pCounter,runCounter,i),4)==66
                        yTest(pCounter,runCounter,i) = 7;
                    elseif info(testingTrials(pCounter,runCounter,i),1)==132 && info(testingTrials(pCounter,runCounter,i),4)==66
                        yTest(pCounter,runCounter,i) = 8;
                    end
                end
            end
            
            
            % Go through each desired time step, as specified at top of program
            for ts = 1:numTimeSteps
                % Build training samples matrix
                % Loop through each frequency band in the tuning data matrix (*.sigChans)
                for fb = 1:numBands
                    % TO SET: Skip 60 Hz band
                    if fb == 10
                        continue;
                    end
                    % If current channel is among those considered "tuned", then
                    % use trials as training samples
                    if max(sigChans{fb,find(tunChans.timeIdxs == timeIdxs(ts),1)}==channels(chan)) > 0
                        X(pCounter,runCounter,:,trainInd(pCounter,runCounter)) = s(freqBands(fb,timeIdxs(ts),trainingTrials(pCounter,runCounter,:)));
                        trainInd(pCounter,runCounter)=trainInd(pCounter,runCounter)+1;
                    end
                end
                % Build test samples matrix
                for fb = 1:numBands
                    % TO SET: Skip 60 Hz band
                    if fb == 10
                        continue;
                    end
                    if max(sigChans{fb,find(tunChans.timeIdxs == timeIdxs(ts),1)}==channels(chan)) > 0
                        Xtest(pCounter,runCounter,:,testInd(pCounter,runCounter)) = s(freqBands(fb,timeIdxs(ts),testingTrials(pCounter,runCounter,:)));
                        testInd(pCounter,runCounter)=testInd(pCounter,runCounter)+1;
                    end
                end
            end
        end
        
        fprintf('Completed p-value %d of %d for channel %d of %d!\n',pCounter,numPs,chan,numChannels);
    end
end

fprintf('Classifying!\n');

% For each p-value and each run
for i = 1:numPs
    for j = 1:numRuns
        % Train Naive Bayes model with X samples matrix
        NBModel =  class_train(s(X(i,j,:,1:trainInd(i,j)-1)),s(y(i,j,:)));
        % Classify samples given in Xtest according to trained model
        yOut(i,j,:) =  class_test(NBModel,s(Xtest(i,j,:,1:testInd(i,j)-1)));
        % Check how we did
        answerCheck(i,j,:,1)=yTest(i,j,:);
        answerCheck(i,j,:,2)=yOut(i,j,:);
        answerCheck(i,j,:,3)=(yTest(i,j,:)==yOut(i,j,:));
        numCorr(i,j)=sum(answerCheck(i,j,:,3));
        percentCorr(i,j)=numCorr(i,j)/numTestTrials;
    end
end

percentChance = 1/length(NBModel.ClassLevels);

% %%%%%%%%%%%% TO SET %%%%%%%%%%%
handle = sprintf('../analysis/workspaces/%s/dirBothClass_50r_160116.mat',folder);
save(handle,'answerCheck','numCorr','percentCorr','percentChance', ...
    'numTrainTrials','numTestTrials','numRuns','numPs','pVals','channels','type');

% Some analysis
meanCorr = mean(percentCorr,2);
devCorr = std(percentCorr,0,2);
figure;
%plot(pVals,meanCorr);
errorbar(pVals,meanCorr,devCorr);
xlabel('p-value');ylabel('Average Prediction %');
title('Pur Direction Prediction Rate vs Confidence Level (50 Runs)');
xlim([.0001 .99]);
ax = gca; set(ax,'XScale','log');
pValues = [.0005 .001, .01, .05, .1, .2, .3, .5, .99];
ax.XTick = pValues;

