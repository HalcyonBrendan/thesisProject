%%Brendan Thorn - 07/15
% anova_code4.m
%
% 27/01/16 - This is a refactored version of anova_code programs that is
% set up for leave-one-out cross-val. IE. It loops through all trials of
% given type, leaving out the i-th one each time, then saves the results
% to hard disk for later use.
%
% This is the same program as anova_code2.m, but accepts input args so it
% can be called from other programs
%
% Runs (one-way) anova on several frequency bands of LFP spectrograms
% to determine if directional tuning exists.
%
% This is a refactored version of anova_code.m.

fprintf('Beginning anova_code4.m!\n');

folder = 'M20100408_256';
dirType = 'trans';
chanType = 'both';%'mip','pmd','both'
% Set which time steps to look at according to *LookSpects matrix
timeIdxs = [61 75 81 88];
numTimeSteps = length(timeIdxs);
% Set a default p-value (this isn't important if you save lots of vars)
pVal = 0.05;

% Load timing and an LFP data structure into workspace for trial info
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

group = info(:,4);
eyeDir = zeros(length(group),1);
eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==11) = 1;
eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==-9) = 2;
eyeDir(Events.EyeHandCoord(:,1)==-10 & Events.EyeHandCoord(:,2)==1) = 3;
eyeDir(Events.EyeHandCoord(:,1)==10 & Events.EyeHandCoord(:,2)==1) = 4;

% Find trials within valid trials for each type of eye movement
validTrials = 1:length(info);
purs = validTrials(group(validTrials)==66 & ~isnan(TT.T66on(validTrials)));
sacs = validTrials(group(validTrials)==65 & ~isnan(TT.T65sac(validTrials)));
cisSacs = sacs(info(sacs,1)>128); transSacs = sacs(info(sacs,1)<129);
cisPurs = purs(info(purs,1)>128); transPurs = purs(info(purs,1)<129);
fixs = validTrials(group(validTrials)==67 & ~isnan(TT.T8(validTrials)));
% saccades which are farther than median from reach
sacThresh = median((TT.T8(sacs)-TT.T65sac(sacs)));
isoSacs = (TT.T8-TT.T65sac) > sacThresh;

numSacs = length(sacs);
numPurs = length(purs);
numCisSacs = length(cisSacs); numTransSacs = length(transSacs);
numCisPurs = length(cisPurs); numTransPurs = length(transPurs);
if strcmp(dirType,'cis')
    dtSacs = cisSacs; dtPurs = cisPurs;
    numDtSacs = numCisSacs; numDtPurs = numCisPurs;
elseif strcmp(dirType,'trans')
    dtSacs = transSacs; dtPurs = transPurs;
    numDtSacs = numTransSacs; numDtPurs = numTransPurs;
end

mipChans = 1:32; pmdChans = 33:48; allChans = 1:48;
if strcmp(chanType,'mip')
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    load(handle);
    channels = mipChans;
elseif strcmp(chanType,'pmd')
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    load(handle);
    chanels = pmdChans;
else % Load relevant channels (all for now)
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    mip=load(handle);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    t = mip.t; f = mip.f;
    clear mip pmd;
    channels = allChans;
end
numChannels = length(channels);

% Set (approximate) bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15; 
bounds(3) = 25; bounds(4) = 35; 
bounds(5) = 45; bounds(6) = 55; 
bounds(7) = 65; bounds(8) = 75; 
bounds(9) = 85; bounds(10) = 95; 
bounds(11) = 105;bounds(12) = 125;
bounds(13) = 150;
numFreqBands = length(bounds);

sacTuningPVals = zeros(numChannels,numDtSacs,numFreqBands,numTimeSteps);
purTuningPVals = zeros(numChannels,numDtPurs,numFreqBands,numTimeSteps);

% Perform anova for the saccades
for chan = 1:numChannels
    % Get correct set of spects for current channel
    spect = s(spects(chan,:,:,:));

    % Average over desired frequency bands
    [freqBands, fBounds] = create_freqBands(spect, f, bounds);
    
    % Now loop through trials. When trial = i, we leave out the i-th trial
    % so it can be used for cross-val later. Sacs first.
    for trial = 1:numDtSacs
        validSacs = [dtSacs(1:trial-1) dtSacs(trial+1:end)];
        % Check for tuning over specified time steps (near eye movement onset)
        for k = 1:numFreqBands
            for j = 1:numTimeSteps
                sacTuningPVals(chan,trial,k,j) = anovan(s(freqBands(k,timeIdxs(j),validSacs)),eyeDir(validSacs,1),'display','off');
            end
        end
    end
    % Purs next.
    for trial = 1:numDtPurs
        validPurs = [dtPurs(1:trial-1) dtPurs(trial+1:end)];
        % Check for tuning over specified time steps (near eye movement onset)
        for k = 1:numFreqBands
            for j = 1:numTimeSteps
                purTuningPVals(chan,trial,k,j) = anovan(s(freqBands(k,timeIdxs(j),validPurs)),eyeDir(validPurs,1),'display','off');
            end
        end
    end
    fprintf('Completed LOO ANOVA for channel %d.\n',channels(chan));
end
% At this point, we have the 4D sacTuningPVals matrix, which tells us for 
% each channel, and for each trial left out for each channel, the p-value 
% at which the frequency bands are directionally tuned as a function of 
% time. Use this matrix to make to create two 3D matrices that 1) give the
% number of tuned channels at each time and frequency band, for each left
% and trial; and 2) list exactly which channels are tuned and comprise the
% number given in the previous matrix. Compare to arbitrarily chosen pVal.
sacSigPVals = s(sum(sacTuningPVals<pVal,1));
sacSigChans = cell(numDtSacs,numFreqBands,numTimeSteps);
for trial = 1:numDtSacs
    for k = 1:numFreqBands
        for j = 1:numTimeSteps
            sacSigChans{trial,k,j}=find(sacTuningPVals(:,trial,k,j)<pVal);
        end
    end
end

purSigPVals = s(sum(purTuningPVals<pVal,1));
purSigChans = cell(numDtPurs,numFreqBands,numTimeSteps);
for trial = 1:numDtPurs
    for k = 1:numFreqBands
        for j = 1:numTimeSteps
            purSigChans{trial,k,j}=find(purTuningPVals(:,trial,k,j)<pVal);
        end
    end
end

times=t(timeIdxs);

handle = sprintf('../analysis/workspaces/%s/%s_dirANOVA_240416.mat',folder,dirType);
save(handle,'dirType','chanType','cisSacs','cisPurs','transSacs','transPurs','sacTuningPVals','purTuningPVals','sacSigChans','purSigChans','sacSigPVals','purSigPVals','times','timeIdxs','fBounds','pVal','sacs','purs','group','eyeDir','folder');

fprintf('Completed anova_code4.m!\n');

