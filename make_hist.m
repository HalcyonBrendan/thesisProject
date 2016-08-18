% make_hist.m  - 12/15
% Brendan Thorn
%
% This code is designed to produce histograms of power at given time in LFP
% spectrograms. Want to check if Naive Bayes Classifier assumption of
% Gaussian distribution or powers is reasonable.

clearvars -except TT

% Specify data set
folder = 'M20100302_645';

% Load an LFP data structure into workspace for trial info
filename = sprintf('../data/monkey_data/LFPChan01.mat');
load(filename);

% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

group = info(:,4);
direction = info(:,1)-128;

% Choose channel of interest
channel = 3;

% Set which time steps to look at
timeIdxs = 271:2:309;
numTimeSteps = length(timeIdxs);

% Set bounds for frequency bands to look at
bounds(1) = 6; bounds(2) = 12; 
bounds(3) = 18; bounds(4) = 24; 
bounds(5) = 30; bounds(6) = 36; 
bounds(7) = 42; bounds(8) = 48; 
bounds(9) = 54; bounds(10) = 66; 
bounds(11) = 80;bounds(12) = 100;
bounds(13) = 150;
numFreqBands = length(bounds);

% Get corresponding spectrogram and t,f data
handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,channel,47);
% Get spectrogram/coherogram/cross-power for channel
[spect,t,f] = getCohs(1,handle,channel,0);
% Center time vector
t=t-7.5;
fprintf('Spectrogram retrieved for channel %d.\n',channel);

% Normalize the spectrogram
%normSpect = normalize_spect(log(permute(spect,[2 1 3])),idxs);
% Or just pretend you normalized
normSpect = log(permute(spect,[2 1 3]));

% Separate normalized spectrograms into bands between 0-200 Hz. Avoid noisy areas ie.
% ~60 Hz and harmonics
[freqBands, fBounds] = create_freqBands(normSpect, f, bounds);

% Find trials for each eye movement type
purs = group==66;
fixs = group==67;
sacs = group==65;
% saccades which are farther than median from reach
sacThresh = median((TT.T8(sacs)-TT.T65sac(sacs)));
isoSacs = (TT.T8-TT.T65sac) > sacThresh;

ups=info(:,1)==129;
downs=info(:,1)==130;
lefts=info(:,1)==131;
rights=info(:,1)==132;

type = sacs&lefts;

% Make histograms for specified trials at a specified time
figure;
xbins1 = -13:.25:-10;
xbins2 = -17:.25:-14;

subplot(2,2,1);
hist(s(freqBands(4,289,type)),xbins1)
xlim([-13 -10])
subplot(2,2,3);
hist(s(freqBands(12,289,type)),xbins2)
xlim([-17 -14])

type = sacs&rights;

subplot(2,2,2);
hist(s(freqBands(4,289,type)),xbins1)
xlim([-13 -10])
subplot(2,2,4);
hist(s(freqBands(12,289,type)),xbins2)
xlim([-17 -14])

sigLevel = [0.95 0.99 0.995];
sigLevel = [1 2 3];
purCorr = [.394 .514 .618];
sacCorr = [.382 .417 .600];


plot(sigLevel, purCorr, sigLevel, sacCorr);
