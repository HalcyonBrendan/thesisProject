%%Brendan Thorn - 11/15
% classification2.m
%
% Just a refactored version of classification.m that allows the use of many
% channels.
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

fprintf('Beginning classification2.m!\n');
clearvars -except TT

numRuns = 10;
confLevels = [5:5:75 77:2:89 90:1:99 991:1:999 9995];
pVals = [.95:-.05:.25 .23:-.02:.11 .10:-.01:.01 .009:-.001:.001 .0005];

numPs = length(pVals);

answerCheck = zeros(100,3,numPs,numRuns);

% Load an LFP data structure into workspace for trial info
%%%%%%%%%% TO SET %%%%%%%%%%%
type = 'sacs'; %purs,both
folder = 'monkey_data';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;
% Get total number of trials from given data
numTrials = size(info,1);

for pCounter = 1:numPs
    
    % Run anova_code2 on first iteration to determine which channels are tuned
    % at a specified confidence level
    if pCounter == 1
        anova_code2_func(pVals(pCounter), confLevels(pCounter));
        %%%%%%%%%%%%%%%% TO SET %%%%%%%%%%%%%%%%%%%%
        % Load tuning information produced by anova_code2.m
        handle = sprintf('../analysis/workspaces/%s/tunedChans_%d_NN_jan0716.mat',folder,confLevels(pCounter));
        tunChans = load(handle);
    else
        for i = 1:size(tunChans.sacTuningPVals,2)
            for j = 1:size(tunChans.sacTuningPVals,3)
                tunChans.sacSigChans{i,j}=find(tunChans.sacTuningPVals(1:numChannels,i,j)<pVals(pCounter));
                tunChans.purSigChans{i,j}=find(tunChans.purTuningPVals(1:numChannels,i,j)<pVals(pCounter));
            end
        end
    end
    
    for runCounter = 1:numRuns   
        %%%%%%%%%%% TODO: Set up proper cross-validation %%%%%%%%%%%
        % Choose approximate proportion of trials you want to train on (difficult
        % to get exact number from here as we will later exclude fixation trials)
        numTrainTrials = round(4/5*numTrials);
        numTestTrials = numTrials-numTrainTrials;
        % Choose which time steps we want to look at for each band
        timeIdxs = 281:2:289;
        numTimeSteps = length(timeIdxs);
        % How many freq bands will we look at
        numBands = length(tunChans.fBounds);
        % Choose which channels we look at at the given time steps
        channels = 1:48;
        numChannels = length(channels);
        
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
        
        %%%%%%%%%%%%%%%%% TO SET %%%%%%%%%%%%%%%%%%%%
        % Randomly select trials that will be used for training
        randTrials = randperm(numTrials);
        trainingTrials = randTrials(1:numTrainTrials);
        % Trials for testing are those not used for training
        testingTrials = randTrials(numTrainTrials+1:end);
        if type == 'sacs'
            trainingTrials = trainingTrials(info(trainingTrials,4)==65);
            testingTrials = testingTrials(info(testingTrials,4)==65);
        elseif type == 'purs'
            trainingTrials = trainingTrials(info(trainingTrials,4)==66);
            testingTrials = testingTrials(info(testingTrials,4)==66);
        else %type == both
            trainingTrials = trainingTrials(info(trainingTrials,4)~=67);
            testingTrials = testingTrials(info(testingTrials,4)~=67);
        end
        % In case repeatability is desired, select first numTrainTrials
        %trainingTrials = 1:numTrainTrials;
        %trainingTrials = trainingTrials(info(trainingTrials,4)==66);
        % In case repeatability is desired, select last numTestTrials
        %testingTrials = numTrainTrials+1:numTrials;
        %testingTrials = testingTrials(info(testingTrials,4)==66);
        
        y = zeros(length(trainingTrials),1);
        yTest = zeros(length(testingTrials),1);
        
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
        
        %%%%%%%%%%%%%%% TO SET %%%%%%%%%%%%%%%%%%%
        % Figure out which channels to use based on which are tuned for both
        % pursuits and saccades, or just choose one or the other
        if type == 'sacs'
            sigChans = tunChans.sacSigChans;
        elseif type == 'purs'
            sigChans = tunChans.purSigChans;
        else
            sigChans = cell(size(tunChans.purSigChans,1),size(tunChans.purSigChans,2));
            for i = 1:size(tunChans.purSigChans,1)
                for j = 1:size(tunChans.purSigChans,2)
                    sigChans{i,j} = intersect(tunChans.purSigChans{i,j},tunChans.sacSigChans{i,j});
                end
            end
        end
        % Simply keep track of index in X sample structures
        trainInd = 1; testInd = 1;
        
        for chan = 1:numChannels
            %%%%%%%%%%%%%%% TO SET %%%%%%%%%%%%%%%%%%
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
            
            % Separate normalized spectrograms into bands between 0-200 Hz. Avoid noisy areas ie.
            % ~60 Hz and harmonics
            [freqBands, fBounds] = create_freqBands(normSpect, f, tunChans.fBounds);
            
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
                        X(:,trainInd) = s(freqBands(fb,timeIdxs(ts),trainingTrials));
                        trainInd=trainInd+1;
                    end
                end
                % Build test samples matrix
                for fb = 1:numBands
                    % TO SET: Skip 60 Hz band
                    if fb == 10
                        continue;
                    end
                    if max(sigChans{fb,find(tunChans.timeIdxs == timeIdxs(ts),1)}==channels(chan)) > 0
                        Xtest(:,testInd) = s(freqBands(fb,timeIdxs(ts),testingTrials));
                        testInd=testInd+1;
                    end
                end
            end
            
            fprintf('Done channel %d of %d\n',chan,length(channels));
        end
        
        fprintf('Training Naive Bayes Classifier!\n');
        
        % Train Naive Bayes model with X samples matrix and their corresponding
        % action types
        NBModel =  class_train(X,y);
        
        fprintf('Classifying test samples!\n');
        
        % Classify samples given in Xtest according to Naive Bayes Model
        % trained above
        yOut =  class_test(NBModel,Xtest);
        
        for idx = 1:length(yTest)
            answerCheck(idx,1,pCounter,runCounter)=yTest(idx);
            answerCheck(idx,2,pCounter,runCounter)=yOut(idx);
            answerCheck(idx,3,pCounter,runCounter)=yTest(idx)==yOut(idx);
        end
        numCorr(pCounter,runCounter)=sum(answerCheck(:,3,pCounter,runCounter));
        percentCorr(pCounter,runCounter)=numCorr(pCounter,runCounter)/length(yOut);
        percentChance(pCounter,runCounter) = 1/length(NBModel.ClassLevels);
        numTrains(pCounter,runCounter) = length(y);
        numTests(pCounter,runCounter) = length(yTest);
        
        clear X Xtest y yTest yOut
        
        fprintf('Completed run %d of %d for p-value %d of %d!\n',runCounter,numRuns,pCounter,numPs);
    end
    if pCounter == 3
        disp('Got here');
    end
end
%%%%%%%%%%%% TO SET %%%%%%%%%%%
save('../analysis/workspaces/monkey_data/dirSacClass_070116.mat','answerCheck','numCorr','percentCorr','percentChance', ...
    'numTrains','numTests','numRuns','numPs','pVals','confLevels');