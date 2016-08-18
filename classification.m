%%Brendan Thorn - 11/15
% classification.m
% 
% Performs classification, determines directional
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

% Get total number of trials from given data
numTrials = size(S1,3);
% Take a guess at how many trials we will train on (note actual number will
% be smaller as we will eliminate fixation trials.
numTrainTrials = 350;
numTestTrials = numTrials-numTrainTrials;
% How many time steps we want to look at for each band
numTimeSteps = 3;
% And how many channels we look at at each time
numChannels = 1;

% Load LFP data structure into workspace for trial info
filename = sprintf('../data/monkey_data/LFPChan01.mat');
load(filename);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Before starting, a bit of code to put eye move times into one vector
alignTimes = zeros(numTrials,1);
for i = 1:numTrials
    if info(i,4)==65
        alignTimes(i)=TT.T65sac(i);
    elseif info(i,4)==66
        alignTimes(i)=TT.T66on(i);
    else
        alignTimes(i)=TT.T8(i);
    end
end

% Separate spectrograms into bands from 0-200 Hz. Avoid noisy areas ie.
% ~60 Hz and harmonics
[freqBands, fBounds] = create_freqBands(permute(S1,[2 1 3]), f, 0);
[freqBands2, fBounds2] = create_freqBands(permute(S2,[2 1 3]), f, 0);

% Normalize freqBands spectrograms to activity around beginning of mem period
% First find correct timing (we need idx of strobe T6)
numNormSteps = 3;
idxTime = (TT.T6-alignTimes)/1000;
idxs = zeros(numTrials,numNormSteps);
for i = 1:numTrials
    idxs(i,1) = find((idxTime(i,:) < t),1);
    idxs(i,2:end) = idxs(i,1)+1:idxs(i,1)+numNormSteps-1;
end

normBands = normalize_spect(freqBands,idxs);
normBands2 = normalize_spect(freqBands2,idxs);

% Concatenate freqBands into a single matrix for ease of use later
if numChannels == 2
    %freqBands = [freqBands; freqBands2];
    normBands = [normBands; normBands2];
end

% Randomly select trials that will be used for training
%randTrials = randperm(numTrials);
randTrials = 1:numTrials;
%trainingTrials = randTrials(1:numTrainTrials);
% In case repeatability is desired, select first numTrainTrials
trainingTrials = 1:numTrainTrials;
trainingTrials = trainingTrials(info(trainingTrials,4)~=67);

% Initialize some structures needed for training
X = zeros(length(trainingTrials),size(normBands,1)*numTimeSteps);
y = zeros(length(trainingTrials),1);

% Build training samples matrix
% X(:,1:3) = squeeze(freqBands(1,291:293,trainingTrials))';
% X(:,4:6) = squeeze(freqBands(2,291:293,trainingTrials))';
% X(:,7:9) = squeeze(freqBands(3,291:293,trainingTrials))';
% X(:,10:12) = squeeze(freqBands(4,291:293,trainingTrials))';
% X(:,13:15) = squeeze(freqBands(5,291:293,trainingTrials))';
% X(:,16:18) = squeeze(freqBands(6,291:293,trainingTrials))';

for i = 1:size(normBands,1)
    X(:,i+(numTimeSteps-1)*(i-1):i+(numTimeSteps-1)*i) = squeeze(normBands(i,275:275+numTimeSteps-1,trainingTrials))';
end


% Build "answers" vector for training samples

% Each combination of up, down, left, right look directions, along with eye
% movement type, saccade or pursuit, will be assigned a numerical value, as
% follows:
%
% Saccade-up     == 1
% Saccade-down   == 2
% Saccade-left   == 3
% Saccade-right  == 4
% Pursuit-up     == 5
% Pursuit-down   == 6
% Pursuit-left   == 7
% Pursuit-right  == 8


for i = 1:length(trainingTrials)
    if info(trainingTrials(i),1)==129 && info(trainingTrials(i),4)==65
        y(i) = 1;
    elseif info(trainingTrials(i),1)==130 && info(trainingTrials(i),4)==65
        y(i) = 2;
    elseif info(trainingTrials(i),1)==131 && info(trainingTrials(i),4)==65
        y(i) = 3;
    elseif info(trainingTrials(i),1)==132 && info(trainingTrials(i),4)==65
        y(i) = 4;
    elseif info(trainingTrials(i),1)==129 && info(trainingTrials(i),4)==66
        y(i) = 5;
    elseif info(trainingTrials(i),1)==130 && info(trainingTrials(i),4)==66
        y(i) = 6;
    elseif info(trainingTrials(i),1)==131 && info(trainingTrials(i),4)==66
        y(i) = 7;
    elseif info(trainingTrials(i),1)==132 && info(trainingTrials(i),4)==66
        y(i) = 8;
    end
end

NBModel =  class_train(X,y);

% Now beging testing

% Trials for testing are those not used for training
testingTrials = randTrials(numTrainTrials+1:end);
testingTrials = testingTrials(info(testingTrials,4)~=67);
% In case repeatability is desired, select last numTestTrials
%testingTrials = numTrainTrials+1:numTrials;

% Initialize some structures needed for testing
Xtest = zeros(length(testingTrials),size(normBands,1)*numTimeSteps);
yTest = zeros(length(testingTrials),1);

% Build testing samples matrix
% Xtest(:,1:3) = squeeze(normBands(1,291:293,testingTrials))';
% Xtest(:,4:6) = squeeze(normBands(2,291:293,testingTrials))';
% Xtest(:,7:9) = squeeze(normBands(3,291:293,testingTrials))';
% Xtest(:,10:12) = squeeze(normBands(4,291:293,testingTrials))';
% Xtest(:,13:15) = squeeze(normBands(5,291:293,testingTrials))';
% Xtest(:,16:18) = squeeze(normBands(6,291:293,testingTrials))';

for i = 1:size(normBands,1)
    Xtest(:,i+(numTimeSteps-1)*(i-1):i+(numTimeSteps-1)*i) = squeeze(normBands(i,275:275+numTimeSteps-1,testingTrials))';
end

% Build "answers" vector for testing samples, as before

for i = 1:length(testingTrials)
    if info(testingTrials(i),1)==129 && info(testingTrials(i),4)==65
        yTest(i) = 1;
    elseif info(testingTrials(i),1)==130 && info(testingTrials(i),4)==65
        yTest(i) = 2;
    elseif info(testingTrials(i),1)==131 && info(testingTrials(i),4)==65
        yTest(i) = 3;
    elseif info(testingTrials(i),1)==132 && info(testingTrials(i),4)==65
        yTest(i) = 4;
    elseif info(testingTrials(i),1)==129 && info(testingTrials(i),4)==66
        yTest(i) = 5;
    elseif info(testingTrials(i),1)==130 && info(testingTrials(i),4)==66
        yTest(i) = 6;
    elseif info(testingTrials(i),1)==131 && info(testingTrials(i),4)==66
        yTest(i) = 7;
    elseif info(testingTrials(i),1)==132 && info(testingTrials(i),4)==66
        yTest(i) = 8;
    end
end

yOut =  class_test(NBModel,Xtest);

answerCheck = zeros(length(yOut),3);
answerCheck(:,1) = yTest; answerCheck(:,2) = yOut; answerCheck(:,3) = yOut==yTest;
numCorr = sum(answerCheck(:,3)); percentCorr = numCorr/length(answerCheck);

