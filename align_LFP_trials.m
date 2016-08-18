%%Brendan Thorn - 11/15
% align_LFP_trials.m
% 
% Aligns LFP time-series trials to desired event (usually strobe in TT.mat
% or equivalent strobe timing file) for delayed eyemove and reach
% experiment.

fprintf('Beginning align_LFP_trials.m!\n');
clear all;

% Load strobe timing file into workspace
folder = 'M20100302_245';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
% Load first LFP channel data into workspace for label info
handle = sprintf('../data/%s/LFPchan01.mat',folder);
load(handle);
info = Events.dataSummary;

% Set some variables to be used later
numChans = 48;
numTrials = size(LFP.AD,1);

% Get label of each trial (look type, then direction)
trialLabel(:,1) = Events.dataSummary(:,4);
trialLabel(:,2) = Events.dataSummary(:,1)-128;

% Create vector of times to be aligned to according to strobes
alignTo = zeros(numTrials,1);
for i = 1:numTrials
    if strcmp(folder,'M20100407_456') && (i==13||i==22||i==163||i==205||i==263||i==273||i==312||i==421)
        alignTo(i)=7500;
    elseif strcmp(folder,'M20100408_256') && (i==179||i==239||i==261||i==267||i==278||i==280)
        alignTo(i)=7500;
    elseif info(i,4)==65
        alignTo(i)=TT.T65sac(i);
    elseif info(i,4)==66
        alignTo(i)=TT.T66on(i);
    else
        alignTo(i)=TT.T8(i);
    end
end
% Instead of aligning to eye movements, align to reach onset
%alignTo = TT.Treach;

% Create matrix to store aligned LFPs
% alignedLFPs = zeros(numTrials,17000,numChans);
% reachAlignedLFPs = zeros(numTrials,17000,numChans);

% For each channel
for i = 1:48
    
    % Load LFP data for channel i
    if i < 10
        handle = sprintf('../data/%s/LFPchan0%d.mat',folder,i);
        load(handle);
    else
        handle = sprintf('../data/%s/LFPchan%d.mat',folder,i);
        load(handle);
    end
    
    % TO SET: put all the trials with NaN alignment times here
    if strcmp(folder,'M20100407_456')
        LFP.badTrials = [13 22 163 205 263 273 312 421];
    elseif strcmp(folder,'M20100408_256')
        LFP.badTrials = [179 239 261 267 278 280];
    end
    
    % Create matrix to store aligned LFPs for channel i
    LFP.alignedLFP = zeros(1,17000);
    %LFP.reachAlignedLFP = zeros(1,17000);
    
    % Align each trial so salient event begins at index 7500
    for j = 1:numTrials
        LFP.alignedLFP(j,7501-alignTo(j):7500-alignTo(j)+length(LFP.AD)) = LFP.AD(j,:);
        %LFP.reachAlignedLFP(j,7501-alignTo(j):7500-alignTo(j)+length(LFP.AD)) = LFP.AD(j,:);
    end
    
    % Add aligned LFP to all-channels matrix
    alignedLFPs(:,:,i) = LFP.alignedLFP;
    %reachAlignedLFPs(:,:,i) = LFP.reachAlignedLFP;
    
    % Save LFP file
    save(handle,'LFP','Events');
    disp(i);
    
end
% 
% % Trim alignedLFPs matrix to minimum size that avoids deleting data
% % Begin by trimming zero-columns before data from any trial begins
% % tempAligned = alignedLFPs;
% tempAligned = reachAlignedLFPs;
% % findEmptyCols = max(squeeze(max(alignedLFPs,[],3)),[],1);
% findEmptyCols = max(squeeze(max(reachAlignedLFPs,[],3)),[],1);
% %for i = 1:length(alignedLFPs)
% for i = 1:length(reachAlignedLFPs)
%     if findEmptyCols(i) > 0
%         tempAligned = tempAligned(:,i:end,:);
%         break;
%     end
% end
% frontColsDeleted = i-1;
% fprintf('Trials aligned at column index %d\n', 7500 - frontColsDeleted);
% % Now trim zero-columns after data from all trials ends
% for i = length(tempAligned):-1:1
%     if findEmptyCols(i) > 0
%         tempAligned = tempAligned(:,1:i,:);
%         break;
%     end
% end
% % alignedLFPs = tempAligned;
% reachAlignedLFPs = tempAligned;
% clear tempAligned;
% % handle = sprintf('../data/%s/alignedLFPs.mat',folder);
% handle = sprintf('../data/%s/reachAlignedLFPs.mat',folder);
% %save(handle,'alignedLFPs','-v7.3');
% %save(handle,'reachAlignedLFPs','-v7.3');
