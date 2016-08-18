%%Brendan Thorn - 11/15
% coherograms.m
% 
% Produces and saves coherograms of LFP data from two given channels. For
% BT M.Eng Thesis, this will primarily be used to compare LFP power signals
% between areas MIP and PMd. Uses cohgramc.m function from Chronux toolbox
% (see chronux.org).
%
% NOTE: Must use LFP data that is properly aligned (especially when working
% across different trials).

fprintf('Beginning coherograms.m!\n');
clear all;

folder = 'M20100302_645';
handle = sprintf('../data/%s/TT.mat',folder);
load(handle);

% Create arrays of channels of interest in MIP and PMd
mipChans = 5;     %1:32;
pmdChans = 36;    %33:48

% Set parameters for cohgramc.m computations later
% Window parameters [length[s] steplength[s]] 
movingwin = [.3 .025];
% Struct with remaining parameters: tapers, pad, Fs, fpass, err, trialave
tapers = [4 7]; pad = 0; Fs = 1000; fpass = [0 Fs/2]; err = [2 1]; trialave = 0;
params = struct('tapers',tapers,'pad',pad,'Fs',Fs,'fpass',fpass,'err',err,'trialave',trialave);

% For channels of interest in MIP (Chans 1-32 in current data (10/15))
for j = 1:length(mipChans)
    % Construct filename of matfile with MIP channel data in it
    if j < 10
        filename = sprintf('../data/%s/LFPChan0%d.mat',folder,mipChans(j));
    else
        filename = sprintf('../data/%s/LFPChan%d.mat',folder,mipChans(j));
    end
    % Load matfile into workspace
    mipChan = load(filename);

    % For channels of interest in PMd (Chans 33-48 in current data (10/15))
    for i = 1:length(pmdChans)

        % Construct filename of matfile with PMd channel data in it
        filename = sprintf('../data/%s/LFPChan%d.mat',folder,pmdChans(i));
        % Load matfile into workspace
        pmdChan = load(filename);

        % Create matrices that stores LFP data to be processed
        data1 = mipChan.LFP.alignedLFP';
        data2 = pmdChan.LFP.alignedLFP';

        [C,phi,S12,S1,S2,t,f,confC,phistd,Cerr] = cohgramc(data1,data2,movingwin,params);

        handle = sprintf('../analysis/coherograms/%s/cohgram_nov1015_chs%d_%d_tpr%d.mat',folder,mipChans(j),pmdChans(i),47);
        save(handle,'C','phi','S12','S1','S2','t','f','confC','phistd','Cerr','-v7.3');
        
        fprintf('MIP Chan %d, PMd Chan %d Complete!\n',mipChans(j),pmdChans(i));
    end
end

% Plot some spectrograms and coherograms

% Bit of code to put alignStrobe times into one vector
for i = 1:length(alignStrobe)
    if mipChan.Events.dataSummary(i,4)==65
        alignTimes(i)=TT.T65sac(i);
    elseif mipChan.Events.dataSummary(i,4)==66
        alignTimes(i)=TT.T66on(i);
    else
        alignTimes(i)=TT.T8(i);
    end
end

% Set look type for plot of interest
type = sacs;
% Set corresponding alignment point
if type == purs
    alignStrobe = TT.T66on;
elseif type == sacs
    alignStrobe = TT.T65sac;
else
    alignStrobe = TT.T8;
end

trialNum = 44;
coh = C(235:355,1:129,trialNum)';
spect1 = S1(235:355,1:129,trialNum)';
spect2 = S2(235:355,1:129,trialNum)';
cSpect = abs(S12(235:355,1:129,trialNum)');

% Display coherograms
figure;
logSpect = log(A1);
% Coherogram (linear scale)
%imagesc(t(1,235:355),f(1,1:129),coh);
% (Cross-)Spectrogram (log scale)
imagesc(t(1,235:355),f(1,1:129),logSpect);
set(gca,'Ydir','Normal');

%draw vertical line at time=0
line([0 0],[0 250]);
%draw vertical line at (average) time of reach
line([(TT.Treach(trialNum)-alignStrobe(trialNum))/1000 (TT.Treach(trialNum)-alignStrobe(trialNum))/1000],[0 250]);
%draw vertical line when (avg) fix point off, signaling reach
%line([(TT.T8(trialNum)-alignStrobe(trialNum))/1000 (TT.T8(trialNum)-alignStrobe(trialNum))/1000],[0 250]);
%draw vertical line when direction cue turns off
line([(TT.T6(trialNum)-alignStrobe(trialNum))/1000 (TT.T6(trialNum)-alignStrobe(trialNum))/1000],[0 250]);
axis([-inf inf 0 250]);

xlabel('Time [s]');ylabel('Frequency [Hz]');
title('Ch.1 Trial 1 DS Spectrogram');
%title('Ch.36 Trial 6 Spectrogram');
%title('Chs.1,36 Trial 11 Cross-Power Spectrogram');
%title('Chs.1,36 Trial 11 Coherogram');
colorbar;

purs=mipChan.Events.dataSummary(:,4)==66;
sacs=mipChan.Events.dataSummary(:,4)==65;
fixs=mipChan.Events.dataSummary(:,4)==67;


% Plot average spectrograms and coherograms
meanSpect1 = mean(S1(215:385,1:129,type),3)';
meanSpect2 = mean(S2(215:385,1:129,type),3)';

% Display plots
figure;
type = purs&rights;
logSpect = log(meanSpectPursRight);
imagesc(t(1,235:355),f(1,1:78),logSpect);
set(gca,'Ydir','Normal');
%draw vertical line at time=0
line([0 0],[0 150]);
%draw vertical line when direction cue turns off
line([mean(TT.T6(type)-alignStrobe(type))/1000 mean(TT.T6(type)-alignStrobe(type))/1000],[0 150]);
axis([-inf inf 0 150]);

xlabel('Time [s]');ylabel('Frequency [Hz]');
title('Ch.2 Mean Right-Pursuit Spectrogram');
caxis([cmin cmax]);
colorbar;

% Plot mean spectrograms for separate look types and directions.
purs=mipChan.Events.dataSummary(:,4)==66;
sacs=mipChan.Events.dataSummary(:,4)==65;
fixs=mipChan.Events.dataSummary(:,4)==67;

ups=mipChan.Events.dataSummary(:,1)==129;
downs=mipChan.Events.dataSummary(:,1)==130;
lefts=mipChan.Events.dataSummary(:,1)==131;
rights=mipChan.Events.dataSummary(:,1)==132;




% Plot average spectrograms and coherograms
meanSpectPursUp = mean(S1(235:355,1:78,purs&ups),3)';
meanSpectPursDown = mean(S1(235:355,1:78,purs&downs),3)';
meanSpectPursLeft = mean(S1(235:355,1:78,purs&lefts),3)';
meanSpectPursRight = mean(S1(235:355,1:78,purs&rights),3)';
meanSpectSacsUp = mean(S1(235:355,1:78,sacs&ups),3)';
meanSpectSacsDown = mean(S1(235:355,1:78,sacs&downs),3)';
meanSpectSacsLeft = mean(S1(235:355,1:78,sacs&lefts),3)';
meanSpectSacsRight = mean(S1(235:355,1:78,sacs&rights),3)';
meanSpectFixsUp = mean(S1(235:355,1:78,fixs&ups),3)';
meanSpectFixsDown = mean(S1(235:355,1:78,fixs&downs),3)';
meanSpectFixsLeft = mean(S1(235:355,1:78,fixs&lefts),3)';
meanSpectFixsRight = mean(S1(235:355,1:78,fixs&rights),3)';

meanSpectPursUp = mean(MS1(1:78,235:355,purs&ups),3);
meanSpectPursDown = mean(MS1(1:78,235:355,purs&downs),3);
meanSpectPursLeft = mean(MS1(1:78,235:355,purs&lefts),3);
meanSpectPursRight = mean(MS1(1:78,235:355,purs&rights),3);
meanSpectSacsUp = mean(MS1(1:78,235:355,sacs&ups),3);
meanSpectSacsDown = mean(MS1(1:78,235:355,sacs&downs),3);
meanSpectSacsLeft = mean(MS1(1:78,235:355,sacs&lefts),3);
meanSpectSacsRight = mean(MS1(1:78,235:355,sacs&rights),3);
meanSpectFixsUp = mean(MS1(1:78,235:355,fixs&ups),3);
meanSpectFixsDown = mean(MS1(1:78,235:355,fixs&downs),3);
meanSpectFixsLeft = mean(MS1(1:78,235:355,fixs&lefts),3);
meanSpectFixsRight = mean(MS1(1:78,235:355,fixs&rights),3);


mmax=max([max(max(log(C1(1:77,235:355)))),max(max(log(C2(1:77,235:355)))),max(max(log(C4(1:77,235:355)))),max(max(log(C6(1:77,235:355)))),max(max(log(C7(1:77,235:355)))),max(max(log(C16(1:77,235:355)))),max(max(log(C25(1:77,235:355)))),max(max(log(C20(1:77,235:355))))]);
mmin=min([min(min(log(C1(1:77,235:355)))),min(min(log(C2(1:77,235:355)))),min(min(log(C4(1:77,235:355)))),min(min(log(C6(1:77,235:355)))),min(min(log(C7(1:77,235:355)))),min(min(log(C16(1:77,235:355)))),min(min(log(C25(1:77,235:355)))),min(min(log(C20(1:77,235:355))))]);
caxis([nmin nmax]);
caxis([cmin cmax]);
caxis([mmin mmax]);