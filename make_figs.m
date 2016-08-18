% Brendan Thorn - 12/15
% make_figs.m - scripts for making plots

%%%%%%%%%%%% Plot of number of tuned channels vs time %%%%%%%%%%%%%%

% Get some colours to use for lines in plot
colors = rand(13,3);
%=[0.0539481928176814,0.548419485731466,0.949370243451058;0.880415938951547,0.0126686341217107,0.647537226254477;0.614411014426690,0.987721046933028,0.0382624597507391;0.00342308456076834,0.0374395360441805,0.123355202520449;0.963402881241132,0.194147564660471,0.293973863482836;0.644568825739417,0.440747029351203,0.680998914337276;0.488862294275121,0.243484263440880,0.644282796413461;0.178710699560875,0.662681633681037,0.187649735135656;0.388925444913291,0.519870901392384,0.416803601643293;0.610385431090712,0.720587737805813,0.782348041459805;0.652039122948022,0.104973057264752,0.112291674642342;0.853291109133284,0.200642858627722,0.413806386178684;0.480308365434127,0.946990212518994,0.776160594864731]
% Get some other values ready
%sacCue = mean(TT.T4(~isnan(TT.T65sac))-TT.T65sac(~isnan(TT.T65sac)));
%purCue = mean(TT.T4(~isnan(TT.T66on))-TT.T66on(~isnan(TT.T66on)));
pVal = 0.01;
sacNumMIPSigPVals = s(sum(sacTuningPVals(1:36,:,:,:)<pVal,1));
sacNumPMdSigPVals = s(sum(sacTuningPVals(37:48,:,:,:)<pVal,1));
purNumMIPSigPVals = s(sum(purTuningPVals(1:36,:,:,:)<pVal,1));
purNumPMdSigPVals = s(sum(purTuningPVals(37:48,:,:,:)<pVal,1));

bands = 1:13; tIdxs = 1:55;
figure;
imagesc(times(tIdxs),bands,s(sacNumMIPSigPVals(1,:,tIdxs))); colorbar('hot');
% for i = 1:length(bands)
%     plot(times(tIdxs), s(sacNumMIPSigPVals(1,bands(i),tIdxs)),'color',colors(i,:));
%     hold on;
% end
% hold off;
% ylim([-.5 25]); xlim([-1 .35]);

% Draw line at direction signal, at last usable time, at eye movement onset
line([0 0], [-.5 36],'LineWidth',1); line([-.149 -.149], [-.5 36],'LineWidth',1); 
%legend('0-5 Hz','5-15 Hz','125-150 Hz','100-125 Hz','35-45 Hz','45-55 Hz');
xlabel('Time Relative to Saccade [s]'); 
%xlabel('Time Relative to Pursuit [s]'); 
ylabel('Num. Tuned Channels /36   (p<0.01)');
title('MIP - ANOVA Saccade Tuning vs Time (p<0.01)');
figure;
for i = 1:length(bands)
    plot(times(tIdxs), s(sacNumPMdSigPVals(1,bands(i),tIdxs)),'color',colors(i,:));
    hold on;
end
hold off;
ylim([-.5 11]); xlim([-1 .35]);
% Draw line at eye movement
line([0 0], [-.5 12],'LineWidth',1); line([-.149 -.149], [-.5 12],'LineWidth',1); 
%legend('0-5 Hz','5-15 Hz','125-150 Hz','100-125 Hz','35-45 Hz','45-55 Hz');
xlabel('Time Relative to Saccade [s]'); 
%xlabel('Time Relative to Pursuit [s]'); 
ylabel('Num. Tuned Channels /12   (p<0.01)');
title('PMd - ANOVA Saccade Tuning vs Time (p<0.01)');
figure;
for i = 1:length(bands)
    plot(times(tIdxs), s(purNumMIPSigPVals(1,bands(i),tIdxs)),'color',colors(i,:));
    hold on;
end
hold off;
ylim([-.5 25]); xlim([-1 .35]);
% Draw line at eye movement
line([0 0], [-.5 36],'LineWidth',1); line([-.149 -.149], [-.5 36],'LineWidth',1); 
%legend('0-5 Hz','5-15 Hz','125-150 Hz','100-125 Hz','35-45 Hz','45-55 Hz');
xlabel('Time Relative to Pursuit [s]'); 
%xlabel('Time Relative to Pursuit [s]'); 
ylabel('Num. Tuned Channels /36   (p<0.01)');
title('MIP - ANOVA Pursuit Tuning vs Time (p<0.01)');
figure;
for i = 1:length(bands)
    plot(times(tIdxs), s(purNumPMdSigPVals(1,bands(i),tIdxs)),'color',colors(i,:));
    hold on;
end
hold off;
ylim([-.5 11]);xlim([-1 .35]);
% Draw line at eye movement
line([0 0], [-.5 12],'LineWidth',1);line([-.149 -.149], [-.5 12],'LineWidth',1);
%legend('0-5 Hz','5-15 Hz','125-150 Hz','100-125 Hz','35-45 Hz','45-55 Hz');
xlabel('Time Relative to Pursuit [s]'); 
%xlabel('Time Relative to Pursuit [s]'); 
ylabel('Num. Tuned Channels /12   (p<0.01)');
title('PMd - ANOVA Pursuit Tuning vs Time (p<0.01)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Figure out how many channels are tuned in each area
pVal = 0.01;
sacNumMIPSigPVals = s(sum(sacTuningPVals(1:36,:,:,:)<pVal,1));
sacNumPMdSigPVals = s(sum(sacTuningPVals(37:48,:,:,:)<pVal,1));

purNumMIPSigPVals = s(sum(purTuningPVals(1:36,:,:,:)<pVal,1));
purNumPMdSigPVals = s(sum(purTuningPVals(37:48,:,:,:)<pVal,1));

%%
%%%%%%%%%%% Plot correct prediction percent vs p-value %%%%%%%%%%%%%

% For more of these, get saved workspace (something like 010116results.mat)
p90both = [0.224;0.286;0.338;0.266;0.171];
p90bothMean = mean(p90both);
p90bothSD = std(p90both);
p90pur = [0.3929;0.5000;0.5294;0.5897;0.4444];
p90purMean = mean(p90pur);
p90purSD = std(p90pur);
p90sac = [0.241;0.533;0.412;0.515;0.548];
p90sacMean = mean(p90sac);
p90sacSD = std(p90sac);

meanSacPercs = [p90sacMean; p95sacMean; p99sacMean; p995sacMean];
meanPurPercs = [p90purMean; p95purMean; p99purMean; p995purMean];
meanBothPercs = [p90bothMean; p95bothMean; p99bothMean; p995bothMean];
errsSac = [p90sacSD; p95sacSD; p99sacSD; p995sacSD];
errsPur = [p90purSD; p95purSD; p99purSD; p995purSD];
errsBoth = [p90bothSD; p95bothSD; p99bothSD; p995bothSD];
X=[.90;.95;.99;.995];

figure;
plot(X,meanSacPercs,X,meanPurPercs,X,meanBothPercs);
set(gca,'XTick',[0.90 0.95 0.99 0.995]);


figure;
errorbar(X,meanSacPercs,errsSac,'-rx'); hold on;
errorbar(X,meanPurPercs,errsPur,'-bx'); hold on;
errorbar(X,meanBothPercs,errsBoth,'-gx');
set(gca,'XTick',[0.90 0.95 0.99 0.995]);
line([.895 1],[.25 .25],'Color','r','LineStyle','--');
line([.895 1],[.25 .25],'Color','b','LineStyle','-.');
line([.895 1],[.125 .125],'Color','g','LineStyle','--');
xlim([.895 1.0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Produce spectrogram
load('/Users/brendan/Documents/Master/Project/LFPs/analysis/coherograms/M20100302_645/spect_dec0215_ch1_tpr47.mat')
t=t-7.5;
spects = permute(S1,[2 1 3]);
trial = 240;
Aspect = s(spects(1:52,255:335,trial));
%Aspect = s(spects(:,:,trial));
logSpect = log(Aspect);
logSpect(isinf(logSpect)) = 0;
%info1 = Events.dataSummary;
%info2 = Events.header;
figure;
imagesc(1000*t(39:88)+149,f(1:78),spect);
colorbar();ax=gca;ax.YDir='normal'; caxis([-16 -10]);
%line([(TT.T4(trial)-TT.T65sac(trial))/1000 (TT.T4(trial)-TT.T65sac(trial))/1000],[f(1) f(78)],'color','c','LineWidth',2);
line([(TT.T6(trial)-TT.T65sac(trial))/1000 (TT.T6(trial)-TT.T65sac(trial))/1000],[f(1) f(78)],'color','c','LineWidth',2);
line([0 0],[f(1) f(78)],'color','g','LineWidth',2);
%line([(TT.T65(trial)-TT.T65sac(trial))/1000 (TT.T65(trial)-TT.T65sac(trial))/1000],[f(1) f(78)],'color','y','LineWidth',2);
line([(TT.Treach(trial)-TT.T65sac(trial))/1000 (TT.Treach(trial)-TT.T65sac(trial))/1000],[f(1) f(78)],'color','r','LineWidth',2);
%line([(TT.T8(trial)-TT.T65sac(trial))/1000 (TT.T8(trial)-TT.T65sac(trial))/1000],[f(1) f(78)],'color','c','LineWidth',2);
%line([(TT.T9(trial)-TT.T65sac(trial))/1000 (TT.T9(trial)-TT.T65sac(trial))/1000],[f(1) f(78)],'color','y','LineWidth',2);
title(trial);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Make histogram of log powers near movement time
folder = 'M20100302_645';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
handle = sprintf('../data/%s/mipLookSpects.mat',folder);
load(handle);
handle = sprintf('../analysis/workspaces/%s/cis_dirANOVA_240416.mat',folder);
tunData = load(handle);

info = Events.dataSummary;

% Select trials to use based on look type
type = 'sacs';
validTrials = 1:length(info);
if strcmp(type,'sacs')
    validTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'purs')
    validTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'both')
    validTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))) & info(validTrials,1)>128 & info(validTrials,1)<133);
end
numTrials = length(validTrials);

% Select channel and freqBand
channel = 13; band = 5;
timeIdxs = [71 73 75 81]';
% Get appropriate spect and average over bands
spect = s(spects(channel,:,:,:));
[freqBands, fBounds] = create_freqBands(spect, f, tunData.fBounds);

upPows = s(freqBands(band,timeIdxs(4),validTrials(info(validTrials,1)==129)));
downPows = s(freqBands(band,timeIdxs(4),validTrials(info(validTrials,1)==130)));
leftPows = s(freqBands(band,timeIdxs(4),validTrials(info(validTrials,1)==131)));
rightPows = s(freqBands(band,timeIdxs(4),validTrials(info(validTrials,1)==132)));
allPows = s(freqBands(band,timeIdxs(4),validTrials));

figure;
subplot(2,2,1);
hist(upPows,8);
title('Up Saccades');
xlabel('Log Power'); ylabel('Count');
subplot(2,2,2);
hist(downPows,8);
title('Down Saccades');
xlabel('Log Power'); ylabel('Count');
subplot(2,2,3);
hist(leftPows,8);
title('Left Saccades');
xlabel('Log Power'); ylabel('Count');
subplot(2,2,4);
hist(rightPows,8);
title('Right Saccades');
xlabel('Log Power'); ylabel('Count');

title('Channel 13, 35-45 Hz Power Histogram');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Make image figure for tuned times and freq bands at two p-vals
folder = 'M20100302_645';
handle = sprintf('../analysis/workspaces/%s/LOO_ANOVA_NN_160216.mat',folder);
tunData = load(handle);

% Select channel and freqBand and pVals
channel = 1; trial = 2;
timeIdxs = [71 73 75]';
pVals = [.05 .01];
numPs = length(pVals);
numBands = length(tunData.fBounds);

A = s(tunData.sacTuningPVals(channel,trial,:,:));
lessTuned = A.*(A<pVals(1));
lessTuned(lessTuned==0) = 1;
moreTuned = A.*(A<pVals(2));
moreTuned(moreTuned==0) = 1;

figure; subplot(2,1,1);
imagesc(tunData.times,1:numBands,lessTuned);c=colorbar;caxis([0 0.06]);
c.Ticks=[.001,.01,.03,.05,.06];c.TickLabels={'p = 0.001','0.01','0.03','0.05','p > T'};
xlabel('Time Relative to Saccade [s]');ylabel('Frequency Bands [Hz]');title('Channel 1, Trial 2 Saccade, Threshold T = 0.05');
ax=gca;ax.YDir='normal'; ax.YTick = 1:13;
ax.YTickLabel = {'0-5','5-15','15-25','25-35','35-45','45-55','55-65','65-75','75-85','85-95','95-105','105-125','125-150'};

subplot(2,1,2);
imagesc(tunData.times,1:numBands,moreTuned);c=colorbar;caxis([0 0.06]);
c.Ticks=[.001,.01,.03,.05,.06];c.TickLabels={'p =0.001','0.01','0.03','0.05','p > T'};
xlabel('Time Relative to Saccade [s]');ylabel('Frequency Bands [Hz]');title('Channel 1, Trial 2 Saccade, Threshold T = 0.01');
ax=gca;ax.YDir='normal'; ax.YTick = 1:13;
ax.YTickLabel = {'0-5','5-15','15-25','25-35','35-45','45-55','55-65','65-75','75-85','85-95','95-105','105-125','125-150'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Compute preferred direction for each channel and plot histogram
type = 'sacs'; %purs,both
folder = 'M20100302_645';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);
handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
handle = sprintf('../analysis/workspaces/%s/LOO_ANOVA_NN_figp99_040316.mat',folder);
tunData = load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Set (approximate) bounds for frequency bands to look at
bounds(1) = 5; bounds(2) = 15; bounds(3) = 25; bounds(4) = 35; bounds(5) = 45; bounds(6) = 55; bounds(7) = 65; 
bounds(8) = 75; bounds(9) = 85; bounds(10) = 95; bounds(11) = 105;bounds(12) = 125;bounds(13) = 150;
numFreqBands = length(bounds);

% Find number of, and particular, appropriate trials based on 'type'
% NOTE: THIS PART MUST SELECT SAME TRIALS AS anova_code4.m !!!
validTrials = 1:length(info);
if strcmp(type,'sacs')
    validTrials = validTrials(info(validTrials,4)==65 & ~isnan(TT.T65sac(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'purs')
    validTrials = validTrials(info(validTrials,4)==66 & ~isnan(TT.T66on(validTrials)) & info(validTrials,1)>128 & info(validTrials,1)<133);
elseif strcmp(type,'both')
    validTrials = validTrials(info(validTrials,4)~=67 & ~(isnan(TT.T65sac(validTrials))&isnan(TT.T66on(validTrials))) & info(validTrials,1)>128 & info(validTrials,1)<133);
end
numTrials = length(validTrials);

chanType = 'both';
channels = 1:48; numChans = length(channels);
if strcmp(chanType,'mip')
    %handle = sprintf('../data/%s/mipLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    load(handle);
elseif strcmp(chanType,'pmd')
    %handle = sprintf('../data/%s/pmdLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    load(handle);
else % Load relevant channels (all for now)
    %handle = sprintf('../data/%s/mipLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/mipLookSpects.mat',folder);
    mip=load(handle);
    %handle = sprintf('../data/%s/pmdLookSpects_100216.mat',folder);
    handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
    pmd=load(handle);
    spects = [mip.spects; pmd.spects];
    t = mip.t; f = mip.f;
    clear mip pmd;
end

% Set time index at which to measure preferred direction (corresponding 
% to tunData.times)
prefIdx = 35;
% And the index corresponding to t in the spectrogram
tIdx = 75;

% Get trials for looks in each direction
us = validTrials(info(validTrials,1)==129); numUs = length(us);
ds = validTrials(info(validTrials,1)==130); numDs = length(ds);
ls = validTrials(info(validTrials,1)==131); numLs = length(ls);
rs = validTrials(info(validTrials,1)==132); numRs = length(rs);

% Get appropriate tuning data
if strcmp(type,'sacs')
    sigChans = s(tunData.sacSigChans(1,:,:));
elseif strcmp(type,'purs')
    sigChans = s(tunData.purSigChans(1,:,:));
end

% Keep track of mean, stddev powers for each direction for each channel
chanPows = zeros(numChans,numFreqBands,4);
chanStds = zeros(numChans,numFreqBands,4);
% And keep track of how many channel-bands are tuned overall, and pref dir
tunedSites = 0; prefDirs = zeros(4,1); minDirs = zeros(4,1);
% Loop through channels
for chan = 1:numChans
    
    % Get log-spectrogram info for correct channel
    spect = s(spects(channels(chan),:,:,:));

    % Average over desired frequency bands
    [freqBands, fBounds] = create_freqBands(spect, f, bounds);
    % For each channel get the mean power response to each direction
    for fb = 1:numFreqBands
        % Check if current channel is tuned in current band
        if max(sigChans{fb,prefIdx}==channels(chan) > 0)
            % Find mean power response to each direction
            chanPows(chan,fb,1) = mean(freqBands(fb,tIdx,us));
            chanPows(chan,fb,2) = mean(freqBands(fb,tIdx,ds));
            chanPows(chan,fb,3) = mean(freqBands(fb,tIdx,ls));
            chanPows(chan,fb,4) = mean(freqBands(fb,tIdx,rs));
            chanStds(chan,fb,1) = std(freqBands(fb,tIdx,us))/sqrt(numUs);
            chanStds(chan,fb,2) = std(freqBands(fb,tIdx,ds))/sqrt(numDs);
            chanStds(chan,fb,3) = std(freqBands(fb,tIdx,ls))/sqrt(numLs);
            chanStds(chan,fb,4) = std(freqBands(fb,tIdx,rs))/sqrt(numRs);
            % Make note that this channel-band was tuned
            tunedSites = tunedSites + 1;
            [~,prefDir] = max(chanPows(chan,fb,:));
            [~,minDir] = min(chanPows(chan,fb,:));
            prefDirs(prefDir) = prefDirs(prefDir) + 1;
            minDirs(minDir) = minDirs(minDir) + 1;
        end
    end

end

% Unroll chanPows, chanStds matrices so they can be plotted
meanPows1 = reshape(chanPows(:,:,1),[],1); stdPows1 = reshape(chanStds(:,:,1),[],1);
meanPows2 = reshape(chanPows(:,:,2),[],1); stdPows2 = reshape(chanStds(:,:,1),[],2);
meanPows3 = reshape(chanPows(:,:,3),[],1); stdPows3 = reshape(chanStds(:,:,1),[],3);
meanPows4 = reshape(chanPows(:,:,4),[],1); stdPows4 = reshape(chanStds(:,:,1),[],4);

meanPows1 = meanPows1(meanPows1~=0); stdPows1 = stdPows1(stdPows1~=0);
meanPows2 = meanPows2(meanPows2~=0); stdPows2 = stdPows2(stdPows2~=0);
meanPows3 = meanPows3(meanPows3~=0); stdPows3 = stdPows3(stdPows3~=0);
meanPows4 = meanPows4(meanPows4~=0); stdPows4 = stdPows4(stdPows4~=0);

figure;
errorbar(meanPows1,stdPows1,'x'); hold on;
errorbar(meanPows2,stdPows2,'o'); hold on;
errorbar(meanPows3,stdPows3,'v'); hold on;
errorbar(meanPows4,stdPows4,'*');


fprintf('Done.\n');

%%
% 3d scatter plots (currently need trainX matrix from svm_LookDetect3)

scatter3(trainX(trainY==0,1),trainX(trainY==0,2),trainX(trainY==0,3));hold on;
scatter3(trainX(trainY==1,1),trainX(trainY==1,2),trainX(trainY==1,3));

