% Make sample LFP-time plot
fig=figure;
plot(-900:325,LFP.alignedLFP(51,6600:7825)/(2*max(abs(LFP.alignedLFP(51,6600:7825)))));
xlim([-900,325]);ylim([-.6 .601]);
line([0 0],[-.6 .6],'color','g','LineWidth',2);
fig.Position = [150 150 570 210];

pVal = 0.05;
numChannels=48;
for i = 1:length(m1.fBounds)
    for j = 1:length(m1.timeIdxs)
        m1.sacSigChans{i,j}=find(m1.sacTuningPVals(1:numChannels,i,j)<pVal);
        m1.purSigChans{i,j}=find(m1.purTuningPVals(1:numChannels,i,j)<pVal);
        m2.sacSigChans{i,j}=find(m2.sacTuningPVals(1:numChannels,i,j)<pVal);
        m2.purSigChans{i,j}=find(m2.purTuningPVals(1:numChannels,i,j)<pVal);
    end
end


line([10^-5 .99], [.25 .25],'LineStyle','--');
xlim([10^-5 .99]);
ylim([.10 .625]);
%%
% Build giant (but not too giant) matrix of relevant parts of the
% spectrograms for all indicated channels so it is unnecessary to
% continually read from disk
folder = 'M20100302_245';
channels = 1:32;
numChans = length(channels);

for chan = 1:numChans
    % Load spectrogram for channels(chan)
    %handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,channels(chan),47);
    handle = '0';
    % Get spectrogram/coherogram/cross-power for channel
    [spect,t,f] = getCohs(1,handle,channels(chan),1,folder);
    % Center time vector
    t=t-7.5; % For normal spect [300ms 25ms]
    %t=t-7.489; % For [240ms 24ms] spect
    fprintf('Spectrogram retrieved for channel %d.\n',channels(chan));
    % Rearrange dimensions and take log
    spect = log(permute(spect,[2 1 3]));
    % Make matrix of interest
    if chan == 1
        % For spect_dec0215 type spects
        spects = zeros(numChans,155,140,size(spect,3));
        % For spect_mar1716 type spects
        %spects = zeros(numChans,93,126,size(spect,3));
    end
    % Now populate it
    spects(chan,:,:,:) = spect(1:155,215:354,:);
    
    fprintf('Done channel %d!\n',channels(chan));
end
% Get appropriate time and freq elements
t=t(215:354);
f=f(1:155);
handle = sprintf('../data/%s/mipLookSpects.mat',folder);
save(handle,'spects','t','f','-v7.3');

clear all;
folder = 'M20100302_245';
channels = 33:48;
numChans = length(channels);
for chan = 1:numChans
    % Load spectrogram for channels(chan)
    %handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,channels(chan),47);
    handle = '0';
    % Get spectrogram/coherogram/cross-power for channel
    [spect,t,f] = getCohs(1,handle,channels(chan),1,folder);
    % Center time vector
    t=t-7.5; % For [300ms 25ms] spect
    %t=t-7.489; % For [240ms 24ms] spect
    fprintf('Spectrogram retrieved for channel %d.\n',channels(chan));
    % Rearrange dimensions and take log
    spect = log(permute(spect,[2 1 3]));
    % Make matrix of interest
    if chan == 1
        % For spect_dec0215 type spects
        spects = zeros(numChans,155,140,size(spect,3));
        % For spect_mar1716 type spects
        %spects = zeros(numChans,93,126,size(spect,3));
    end
    % Now populate it
    spects(chan,:,:,:) = spect(1:155,215:354,:);
    
    fprintf('Done channel %d!\n',channels(chan));
end
% Get appropriate time and freq elements
t=t(215:354);
f=f(1:155);
handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
save(handle,'spects','t','f','-v7.3');
% %%
% % Forgot to save t and f vectors with spects, so reopen, add, resave
% folder = 'M20100407_456';
% 
% handle = sprintf('../data/%s/mipLookSpects.mat',folder);
% load(handle);
% 
% % Load spectrogram for channels(1) just to see t and f
% handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,37,47);
% % Get spectrogram/coherogram/cross-power for channel
% [spect,t,f] = getCohs(1,handle,37,0);
% % Center time vector
% t=t-7.5;
% 
% % Modify t and f vectors to match spects
% time=t;
% freq=f;
% t = t(201:400);
% f = f(1:93);
% 
% handle = sprintf('../data/%s/pmdLookSpects.mat',folder);
% save(handle,'spects','t','f','-v7.3');
% %%
% oneCorr=sum(check(:,1)==1 & check(:,3)==1)/sum(check(:,1)==1);
% twoCorr=sum(check(:,1)==2 & check(:,3)==1)/sum(check(:,1)==2);
% threeCorr=sum(check(:,1)==3 & check(:,3)==1)/sum(check(:,1)==3);
% fourCorr=sum(check(:,1)==4 & check(:,3)==1)/sum(check(:,1)==4);
% 
% 
% %%
% % 030116 BT
% %
% % Check total power on trials for normalization purposes
% spect = s(spects(1,:,:,:));
% 
% for i = 100:105
%    figure;
%    plot(t,sum(s(spect(:,:,valTrials(i))),1));
%    xlim([-3 3]);ylim([-1400 -1100]);
% end
% %%
% % 040116 BT
% %
% % Check whether total power varies systematically with look direction
% for i = 1:36
%     spect = s(spects(i,:,:,:));
% 
%     dirPow(i,1) = mean(sum(spect(:,89,valTrials(info(valTrials,1)==129)),1));
%     stdErr(i,1) = std(sum(spect(:,89,valTrials(info(valTrials,1)==129)),1))/sqrt(length(valTrials(info(valTrials,1)==129)));
%     dirPow(i,2) = mean(sum(spect(:,89,valTrials(info(valTrials,1)==130)),1));
%     stdErr(i,2) = std(sum(spect(:,89,valTrials(info(valTrials,1)==130)),1))/sqrt(length(valTrials(info(valTrials,1)==130)));
%     dirPow(i,3) = mean(sum(spect(:,89,valTrials(info(valTrials,1)==131)),1));
%     stdErr(i,3) = std(sum(spect(:,89,valTrials(info(valTrials,1)==131)),1))/sqrt(length(valTrials(info(valTrials,1)==131)));
%     dirPow(i,4) = mean(sum(spect(:,89,valTrials(info(valTrials,1)==132)),1));
%     stdErr(i,4) = std(sum(spect(:,89,valTrials(info(valTrials,1)==132)),1))/sqrt(length(valTrials(info(valTrials,1)==132)));
% end
% 
% % Check whether total power varies systematically with look direction
% 
% spect = s(spects(22,:,:,:));
% 
% sortPow1 = s(sort(sum(spect(:,89,valTrials(info(valTrials,1)==129)),1)));
% sortPow2 = s(sort(sum(spect(:,89,valTrials(info(valTrials,1)==130)),1)));
% sortPow3 = s(sort(sum(spect(:,89,valTrials(info(valTrials,1)==131)),1)));
% sortPow4 = s(sort(sum(spect(:,89,valTrials(info(valTrials,1)==132)),1)));

%%
for i = 1:48
    if i==27
        continue;
    end
spect = s(spects(channels(8),:,:,:));
[freqBands, fBounds] = create_freqBands(spect, f, bounds);
data1 = s(freqBands(5,81,validSacs));
data2 = s(freqBands(5,74,validSacs));
sdev1 = std(data1,0,1);
sdev2 = std(data2,0,1);
vect1 = reshape(data1,[],1);
vect2 = reshape(data2,[],1);

figure;
scatter(1,mean(data1,1),'bo');hold on;
scatter(1,mean(data1,1)+sdev1,'go');hold on;
scatter(1,mean(data1,1)-sdev1,'go');hold on;
scatter(1,mean(data2,1),'rx');hold on;
scatter(1,mean(data2,1)+sdev2,'cx');hold on;
scatter(1,mean(data2,1)-sdev2,'cx');hold off;

figure;
plot(vect1,'bo'); hold on;
plot(vect2,'rx'); hold on;

end

%%
% Smooth cursor movements

currCur = yCur.LFP;
currCur.smooth = zeros(size(currCur.AD));

for i = 1:size(currCur.AD,1)
    currCur.smooth(i,:) = smooth(currCur.AD(i,:)',20)';
end
yCur.LFP = currCur;

%%
% For direction analysis

handle = sprintf('../data/%s/LFPChan01.mat',folder);
load(handle);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Get types of trials
transSacs = find(info(:,4)==65 & info(:,1)<129);
transSacIdxs = find(ismember(validTrials,transSacs))';
transPurs = find(info(:,4)==66 & info(:,1)<129);
transPurIdxs = find(ismember(validTrials,transPurs));
cisSacs = find(info(:,4)==65 & info(:,1)>128);
cisSacIdxs = find(ismember(validTrials,cisSacs))';
cisPurs = find(info(:,4)==66 & info(:,1)>128);
cisPurIdxs = find(ismember(validTrials,cisPurs));

compressedAnswers = s(answerCheck(:,:,:,:,3));

meanTransSacs = s(mean(compressedAnswers(:,:,transSacIdxs),3));
meanCisSacs = s(mean(compressedAnswers(:,:,cisSacIdxs),3));
meanTransPurs = s(mean(compressedAnswers(:,:,transPurIdxs),3));
meanCisPurs = s(mean(compressedAnswers(:,:,transPurIdxs),3));

epMeanTransSacs = mean(meanTransSacs,2);
epMeanCisSacs = mean(meanCisSacs,2);
epMeanTransPurs = mean(meanTransPurs,2);
epMeanCisPurs = mean(meanCisPurs,2);
