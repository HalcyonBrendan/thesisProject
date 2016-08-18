%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Key variables to load before getting started
folder = 'M20100302_645';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
info = Events.dataSummary;
handle = sprintf('../data/%s/mipLookSpects',folder);
load(handle);
% mip=load(handle);
% handle = sprintf('../data/%s/pmdLookSpects',folder);
% pmd=load(handle);
% spects = [mip.spects; pmd.spects];
% t = mip.t;f=mip.f;
% clear mip pmd;
channels = 1:36; numChans = length(channels);

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
% 010216 BT
% Get SSE between power signals on unseen trials and mean power at look
% time on training trials. Use a bunch of channels.
numCompares = 4;
testSig = zeros(numFreqBands,length(t),numCompares);
numTrains = 110; maxTrial = 120;
trainIdx = numTrains-1;
for testTrial = numTrains+1:maxTrial
    trainIdx = trainIdx+1;
    SSEcum = zeros(1,length(t));
    SSE = zeros(1,length(t));
    for chan = 1:numChans
        spect = s(spects(chan,:,:,:));
        [freqBands, fBounds] = create_freqBands(spect, f, bounds);
        %normBands = normByBand(freqBands,1);
        meanPows = mean(freqBands(:,85:88,valTrials(1:trainIdx)),3);
        for c = 1:numCompares
            testSig(:,:,c) = [s(freqBands(:,c:end,valTrials(testTrial))) ones(numFreqBands,c-1) ];
        end
        SE = (testSig-permute(repmat(meanPows,1,1,length(t)),[1 3 2])).^2;
        % Ignore some bands that contribute too much noise at look time
        SE(1:10,:,:)=0;
        SSE = s(sum(sum(SE,1),3));
        SSE(isinf(SSE))=-.1;
        SSEcum = SSEcum+SSE.*(SSEcum>=0);
        SSEcum(SSEcum>(5-3*chan/numChans)*mean(SSEcum(SSEcum>0),2))=-.1;
    end
    SSEcum = SSEcum/numChans;
    [errs,idxs]=sort(SSEcum);
    [minErr,minIdx] = find(errs>0,1); 
    titl = sprintf('Saccade Detection by Error. Trial %d',valTrials(testTrial));
    figure;plot(t,SSEcum);
    xlabel('Time rel. Saccade');ylabel('Power SSE rel. Saccade Mean');title(titl);
    dim = [.2 .5 .3 .3];
    str = sprintf('Minimum Error at t = %f',t(idxs(minIdx)));
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    fprintf('Done trial %d\n',testTrial);
end


%%
% 010216 BT
% Determine which frequency bands are helpful or unhelpful for determining
% look timing
numTrains = 150; maxTrial = 160;
trainIdx = numTrains-1;
SEcum = zeros(maxTrial-numTrains,numFreqBands,length(t));
for testTrial = numTrains+1:maxTrial
    trainIdx = trainIdx+1;
    SE = zeros(numFreqBands,length(t));
    for chan = 1:numChans
        spect = s(spects(chan,:,:,:));
        [freqBands, fBounds] = create_freqBands(spect, f, bounds);
        %normBands = normByBand(freqBands,1);
        medPows = median(freqBands(:,89,valTrials(1:trainIdx)),3);
        testSig = freqBands(:,:,valTrials(testTrial));
        %SSE = sum((testSig-repmat(medPows,1,length(t))).^2,1);
        SE = (testSig-repmat(medPows,1,length(t))).^2;
        SE(isinf(SE))=-.1;
        SEcum(testTrial-numTrains,:,:) = s(SEcum(testTrial-numTrains,:,:))+SE.*(s(SEcum(testTrial-numTrains,:,:))>=0);
        %SEcum(SEcum>(5-3*chan/numChans)*mean(SEcum(SEcum>0)))=-.1;
    end
    SEcum(testTrial-numTrains,:,:) = SEcum(testTrial-numTrains,:,:)/numChans;
    titl = sprintf('Saccade Detection by Error. Trial %d',valTrials(testTrial));
    figure;
    for fb = 1:numFreqBands
        plot(t(25:200),s(SEcum(testTrial-numTrains,fb,25:200)));
        xlabel('Time rel. Saccade');ylabel('Power SSE rel. Saccade Mean');title(titl);
        hold on;
    end
    hold off;
    fprintf('Done trial %d\n',testTrial);
end
%%
% 020216 BT
% Try using SVM to classify eye movement vs not.
fprintf('Using SVM to detect eye movement times.\n');

%Get some details in order
type = 'sacs'; %%%%%%%%%%%%%%%%%%%%%%% <--- SET THIS 
if strcmp(type,'sacs')
    valTrials = validSacTrials;
elseif strcmp(type,'purs')
    valTrials = validPurTrials;
else
    valTrials = validBothTrials;
end
numTrials = length(valTrials);

% Choose which trials to train and test on (idxs within valTrials)
trainTrials = (1:170)';
testTrials = (171)';
numTrains = length(trainTrials);numTests = length(testTrials);
% Choose two times which will be used as "MOVE=1" 
move = [86; 87; 88; 89; 90]; numMoves = length(move);
% Select four others (initially) as "NO MOVE=2"
% If you want, select randomly
noMove = [30:3:78 99:3:180]'; numNoMoves = length(noMove);

% Declare svm matrices
svmTrainingX = zeros(numTrains*(numMoves+numNoMoves),numChans*numFreqBands);
svmTestingX = zeros(numTests*(numMoves+numNoMoves),numChans*numFreqBands);
svmTrainingY = zeros(numTrains*(numMoves+numNoMoves),1);
svmTestingY = zeros(numTests*(numMoves+numNoMoves),1);
trainRowInd = 1; testRowInd = 1; colInd = 1;
% Loop through channels
numChans = size(spects,1);
for chan = 1:numChans
    % reset rowInds
    trainRowInd = 1; testRowInd = 1;
    % Separate spectrogram into frequency bands
    spect = s(spects(chan,:,:,:));
    [freqBands, fBounds] = create_freqBands(spect, f, bounds);
    % Get "MOVE" examples from this channel
    for i = 1:length(move)
        % Get data for specified time for training
        moveX = s(freqBands(:,move(i),valTrials(trainTrials)))';
        % Append it to training matrix
        svmTrainingX(trainRowInd:trainRowInd+size(moveX,1)-1,colInd:colInd+size(moveX,2)-1) = moveX;
        if chan == 1
            svmTrainingY(trainRowInd:trainRowInd+size(moveX,1)-1) = ones(size(moveX,1),1);
        end
        % Update indices
        trainRowInd = trainRowInd + size(moveX,1);
        % Get data for specified time for testing (only one "move" time)
        if i == 3
            moveX = s(freqBands(:,move(i),valTrials(testTrials)))';
            % Append it to testing matrix
            svmTestingX(testRowInd:testRowInd+size(moveX,1)-1,colInd:colInd+size(moveX,2)-1) = moveX;
            if chan == 1
                svmTestingY(testRowInd:testRowInd+size(moveX,1)-1) = ones(size(moveX,1),1);
            end
            % Update indices
            testRowInd = testRowInd + size(moveX,1);
        end
    end
    % Get "NO MOVE" examples from this channel
    for i = 1:length(noMove)
        % Get data for specified time for training
        noMoveX = s(freqBands(:,noMove(i),valTrials(trainTrials)))';
        % Append it to training matrix
        svmTrainingX(trainRowInd:trainRowInd+size(noMoveX,1)-1,colInd:colInd+size(noMoveX,2)-1) = noMoveX;
        if chan == 1
            svmTrainingY(trainRowInd:trainRowInd+size(noMoveX,1)-1) = zeros(size(noMoveX,1),1);
        end
        % Update indices
        trainRowInd = trainRowInd + size(noMoveX,1);
        % Get data for specified time for testing
        noMoveX = s(freqBands(:,noMove(i),valTrials(testTrials)))';
        % Append it to training matrix
        svmTestingX(testRowInd:testRowInd+size(noMoveX,1)-1,colInd:colInd+size(noMoveX,2)-1) = noMoveX;
        if chan == 1
            svmTestingY(testRowInd:testRowInd+size(noMoveX,1)-1) = zeros(size(noMoveX,1),1);
        end
        % Update indices
        testRowInd = testRowInd + size(noMoveX,1);
    end
    colInd = colInd + size(noMoveX,2);
end

% Only use the relevant parts of the matrices
svmTrainingX = svmTrainingX(1:trainRowInd-1,1:colInd-1);
svmTestingX = svmTestingX(1:testRowInd-1,1:colInd-1);
svmTrainingY = svmTrainingY(1:trainRowInd-1);
svmTestingY = svmTestingY(1:testRowInd-1);

% Train classifier
SVMModel = fitcsvm(svmTrainingX,svmTrainingY,'Standardize',true,'Prior',[.99 .01],'BoxConstraint',1);
% Cross-validate
%CVSVMModel = crossval(SVMModel);
% Save some space
%SVMModel = compact(SVMModel);
% Fit posterior function for prediction
SVMModel = fitPosterior(SVMModel,svmTrainingX,svmTrainingY);
% Predict on testing set
[svmTestingYout,score] = predict(SVMModel,svmTestingX);
% Display results
results = zeros(testRowInd-1,4);
results(:,1) = svmTestingY;results(:,2)=svmTestingYout;results(:,3)=(results(:,1)==results(:,2));results(:,4)=score(:,2);
moveCorr=sum(results(results(:,1)==1,3)==1)/sum(results(:,1)==1);noMoveCorr=sum(results(results(:,1)==0,3)==1)/sum(results(:,1)==0);
%results = table(svmTestingY,svmTestingYout,score(:,2),'VariableNames',{'TrueLabels','PredictedLabels','PosClassPosterior'});

fprintf('Completed SVM detection.\n');

%%
% MSE stuff
purMSEf = zeros(numChans,363-235+1);
purMSEs = zeros(numChans,363-235+1);
sacMSEf = zeros(numChans,363-235+1);
sacMSEs = zeros(numChans,363-235+1);
purSEf = zeros(numFreqBands,363-235+1);
purSEs = zeros(257,363-235+1);
sacSEf = zeros(numFreqBands,363-235+1);
sacSEs = zeros(257,363-235+1);

for chan = 1:numChans

    % Plot power in freq bands vs time

    % Load spectrogram
    handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,channels(chan),47);
    [spect,t,f] = getCohs(1,handle,channels(1),channels(1),folder);
    % Center time vector
    t=t-7.5;
    fprintf('Spectrogram retrieved for channel %d.\n',channels(chan));
    
    %     % (OPTIONAL) Normalize the spectrogram
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
    % Or don't normalize, but pretend you did so code still works
    normSpect = log(permute(spect,[2 1 3]));
    % Get rid of infs created by log to avoid problems later
    normSpect(isinf(normSpect)) = 0;
    
    % Average over desired frequency bands
    % Separate normalized spectrograms into bands between 0-200 Hz. Avoid noisy areas ie.
    % ~60 Hz and harmonics
    [freqBands, fBounds] = create_freqBands(normSpect, f, bounds);
    
    trialSacAv=mean(freqBands(:,:,validSacTrials),3);
    trialPurAv=mean(freqBands(:,:,validPurTrials),3);
    
    spectSacAv=mean(normSpect(:,:,validSacTrials),3);
    spectPurAv=mean(normSpect(:,:,validPurTrials),3);
    
    % Compute mean square differences between averaged spects and spects of
    % left out trials vs time (purs 80,81, sacs 82,83 on 25/10/16)
    %purMSEf(chan,:) = mean((trialPurAv(:,235:363)-s(freqBands(:,235:363,81))).^2,1);
    %purMSEs(chan,:) = mean((spectPurAv(:,235:363)-s(normSpect(:,235:363,81))).^2,1);
    %sacMSEf(chan,:) = mean((trialPurAv(:,235:363)-s(freqBands(:,235:363,82))).^2,1);
    %sacMSEs(chan,:) = mean((spectPurAv(:,235:363)-s(normSpect(:,235:363,82))).^2,1);
    %purSEf = purSEf+(trialPurAv(:,235:363)-s(freqBands(:,235:363,81))).^2;
    %purSEs = purSEs+(spectPurAv(:,235:363)-s(normSpect(:,235:363,81))).^2;
    sacSEf = sacSEf+(trialPurAv(:,235:363)-s(freqBands(:,235:363,83))).^2;
    %sacSEs = sacSEs+(spectPurAv(:,235:363)-s(normSpect(:,235:363,82))).^2;
    
    % Plot 
%     figure;
%     plot(t(235:363),purMSE);
%     figure;
%     plot(t(235:363),sacMSE);

%     for fb = 1:length(fBounds)
%         figure;
%         plot(t(235:363),trialSacAv(fb,235:363));
%         xlabel('Time Rel. to Eye Movement');ylabel('Mean Log power in Band');
%         titl = sprintf('Ch. %d Power vs. Time in Freq. Band <%f Hz.',channels(chan),fBounds(fb));
%         title(titl);
%     end
%     
%     for fb = 1:length(fBounds)
%         figure;
%         plot(t(235:363),trialPurAv(fb,235:363));
%         xlabel('Time Rel. to Eye Movement');ylabel('Log power in Band');
%         titl = sprintf('Ch. %d Power vs. Time in Freq. Band <%f Hz.',channels(chan),fBounds(fb));
%         title(titl);
%     end
    
end

for fb = 1:length(fBounds)
%     figure;
%     plot(t(235:363),purSEf(fb,:));
%     xlabel('Time Rel. to Eye Movement');ylabel('SE in Band');
%     titl = sprintf('PMd Pur SE vs. Time in Freq. Band <%f Hz.',fBounds(fb));
%     title(titl);
    figure;
    plot(t(235:363),sacSEf(fb,:));
    xlabel('Time Rel. to Eye Movement');ylabel('MSE');
    if fb > 1
        titl = sprintf('Trial 83 MIP Sac MSE (rel to avg power for sac) vs. Time in Freq. Band %f < f <%f Hz.',fBounds(fb-1),fBounds(fb));
    else
        titl = sprintf('Trial 83 PMd Sac MSE vs. Time in Freq. Band f <%f Hz.',fBounds(fb));
    end
    title(titl);
end

% figure;
% plot(t(235:363),mean(purMSEf,1));
% title('PMd Freq Pur 81 MSE');
% figure;
% plot(t(235:363),mean(purMSEs,1));
% title('PMd Spect Pur 81 MSE');
% figure;
% plot(t(235:363),mean(sacMSEf,1));
% title('PMd Freq Sac 83 MSE');
% figure;
% plot(t(235:363),mean(sacMSEs,1));
% title('PMd Spect Sac 83 MSE');


% Get some info about power in various bands at look time
spect = s(spects(1,:,:,:));
[freqBands, fBounds] = create_freqBands(spect, f, bounds);

trainMeanSig = zeros(numFreqBands,length(t));
testSig1 = zeros(numFreqBands,length(t));
testSig2 = zeros(numFreqBands,length(t));
for fb = 1:numFreqBands
    trainMeanSig(fb,:) = mean(freqBands(fb,:,validSacTrials(1:numTrains)),3);
    testSig1(fb,:) = freqBands(fb,:,validSacTrials(101));
    testSig2(fb,:) = freqBands(fb,:,validSacTrials(102));
end
SSE1 = sum((testSig1-repmat(meanPows,1,length(t))).^2,1);
SSE2 = sum((testSig2-repmat(meanPows,1,length(t))).^2,1);
SSE1(isinf(SSE1))=0;SSE2(isinf(SSE2))=0;
figure;plot(t(30:end),SSE1(30:end));figure;plot(t(30:end),SSE2(30:end));



numTrains = 100;
meanPows = mean(freqBands(:,89,validSacTrials(1:numTrains)),3);
medPows = median(freqBands(:,89,validSacTrials(1:numTrains)),3);
stdPows = std(freqBands(:,89,validSacTrials(1:numTrains)),0,3);
minPows = min(freqBands(:,89,validSacTrials(1:numTrains)),[],3);
maxPows = max(freqBands(:,89,validSacTrials(1:numTrains)),[],3);
snrPows = meanPows./stdPows;
statPows(:,1)=maxPows;
statPows(:,2)=medPows;
statPows(:,3)=minPows;
statPows(:,4)=meanPows;
statPows(:,5)=stdPows;

[S,I] = sort(freqBands(:,89,validSacTrials),3);
S=s(S);I=s(I);

pows = zeros(numFreqBands,length(validSacTrials)-numTrains);
checklist = zeros(numFreqBands,length(validSacTrials)-numTrains);
for i = numTrains+1:length(validSacTrials)
    pows(:,i-numTrains) = freqBands(:,89,validSacTrials(i));
    checklist(:,i-numTrains) = (pows(:,i-numTrains)<statPows(:,2)+statPows(:,5)) & (pows(:,i-numTrains)>statPows(:,2)-statPows(:,5));
end
checksum = sum(checklist,1);

