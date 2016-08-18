%%Brendan Thorn - 07/15
% anova_code.m
% 
% Runs (one-way) anova on several frequency bands of LFP spectrograms
% to determine if directional tuning exists.
%
% Note: Before running, load LFP data (eg. LFPChan01.mat - only need
% one channel to get trial info) and strobe timing data (eg. TT.mat)
% into workspace. Also adjust handles within code as necessary.

memTimes = -70:4:90;
memTimes(2,:) = memTimes(1,:)*.025;
group = Events.dataSummary(:,4);
direction = Events.dataSummary(:,1)-128;

purs = group==66;
fixs = group==67;
sacs = group==65;
% saccades which are farther than median from reach
sacThresh = median((TT.T8(sacs)-TT.T65sac(sacs)));
isoSacs = (TT.T8-TT.T65sac) > sacThresh;

% ****** TO BE SET ****************************************************
% set trials to whichever eye movement type you are interested in *****
trials = purs;

tuningPVals1 = zeros(48,41);
tuningPVals2 = zeros(48,41);
tuningPVals3 = zeros(48,41);
tuningPVals4 = zeros(48,41);
chanArray = cell(48,1);

%for each channel
for i = 1:48
    %every 4th channel, load next four channels
    if mod(i,4)==1
        %handle = sprintf('../analysis/workspaces/M20100301_246/chanArray%d-%dAug0315.mat',i,i+3);
        handle = sprintf('../analysis/workspaces/monkey_data/chanArray%d-%dJul1715.mat',i,i+3);
    end
    clear chanArray;
    load(handle);
    for j = 1:length(memTimes)
        %consider mean power in ~0-10Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(1:5,221-6+memTimes(1,j),:)));
        tuningPVals1(i,j) = anovan(temp(trials),direction(trials,1));
        close(1);
        %consider mean power in ~10-20Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(6:11,221-6+memTimes(1,j),:)));
        tuningPVals2(i,j) = anovan(temp(trials),direction(trials,1));
        close(1);
        %consider mean power in ~20-30Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(12:16,221-6+memTimes(1,j),:)));
        tuningPVals3(i,j) = anovan(temp(trials),direction(trials,1));
        close(1);
        %consider mean power in ~30-40Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(17:21,221-6+memTimes(1,j),:)));
        tuningPVals4(i,j) = anovan(temp(trials),direction(trials,1));
        close(1);
        %consider mean power in ~40-60Hz freq band
        temp = squeeze(mean(chanArray{i,1}.spectsAlign(22:30,221-6+memTimes(1,j),:)));
        tuningPVals5(i,j) = anovan(temp(trials),direction(trials,1));
        close(1);
    end
    disp(i);
end

% %FOR SACCADE TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot figure for MIP channels
% sigPVals1 = sum(tuningPVals1(1:32,:)<.01);sigPVals2 = sum(tuningPVals2(1:32,:)<.01);
% sigPVals3 = sum(tuningPVals3(1:32,:)<.01);sigPVals4 = sum(tuningPVals4(1:32,:)<.01);
% sigPVals5 = sum(tuningPVals5(1:32,:)<.01);
% figure;
% X = memTimes(2,:)*1000;
% plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to saccade [ms]');ylabel('Number of significant channels (p<.01) /32');
% legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');title('MIP Saccades');
% %draw vertical line at time=0
% line([mean(TT.T65sac(trials)-TT.T65sac(trials)) mean(TT.T65sac(trials)-TT.T65sac(trials))],[0 32]);
% %draw vertical line at (average) time of reach
% line([mean(TT.Treach(trials)-TT.T65sac(trials)) mean(TT.Treach(trials)-TT.T65sac(trials))],[0 32]);
% %draw vertical line when (avg) fix point off, signaling reach
% line([mean(TT.T8(trials)-TT.T65sac(trials)) mean(TT.T8(trials)-TT.T65sac(trials))],[0 32]);
% %draw vertical line when direction cue turns on
% line([mean(TT.T4(trials)-TT.T65sac(trials)) mean(TT.T4(trials)-TT.T65sac(trials))],[0 32]);
% axis([-inf inf 0 32]);
% %draw vertical line when direction cue turns off
% line([mean(TT.T6(trials)-TT.T65sac(trials)) mean(TT.T6(trials)-TT.T65sac(trials))],[0 32]);
% axis([-inf inf 0 32]);
% %plot figure for PMd channels
% sigPVals1 = sum(tuningPVals1(33:48,:)<.01);sigPVals2 = sum(tuningPVals2(33:48,:)<.01);
% sigPVals3 = sum(tuningPVals3(33:48,:)<.01);sigPVals4 = sum(tuningPVals4(33:48,:)<.01);
% sigPVals5 = sum(tuningPVals5(33:48,:)<.01);
% figure;
% X = memTimes(2,:)*1000;
% plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to saccade [ms]');ylabel('Number of significant channels (p<.01) /16');
% legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');title('PMd Saccades');
% %draw vertical line at time=0
% line([mean(TT.T65sac(trials)-TT.T65sac(trials)) mean(TT.T65sac(trials)-TT.T65sac(trials))],[0 16]);
% %draw vertical line at (average) time of reach
% line([mean(TT.Treach(trials)-TT.T65sac(trials)) mean(TT.Treach(trials)-TT.T65sac(trials))],[0 16]);
% %draw vertical line when (avg) fix point off, signaling reach
% line([mean(TT.T8(trials)-TT.T65sac(trials)) mean(TT.T8(trials)-TT.T65sac(trials))],[0 16]);
% %draw vertical line when direction cue turns on
% line([mean(TT.T4(trials)-TT.T65sac(trials)) mean(TT.T4(trials)-TT.T65sac(trials))],[0 16]);
% axis([-inf inf 0 16]);
% %draw vertical line when direction cue turns off
% line([mean(TT.T6(trials)-TT.T65sac(trials)) mean(TT.T6(trials)-TT.T65sac(trials))],[0 16]);
% axis([-inf inf 0 16]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FOR PURSUIT TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot figure for MIP channels
sigPVals1 = sum(tuningPVals1(1:32,:)<.01);sigPVals2 = sum(tuningPVals2(1:32,:)<.01);
sigPVals3 = sum(tuningPVals3(1:32,:)<.01);sigPVals4 = sum(tuningPVals4(1:32,:)<.01);
sigPVals5 = sum(tuningPVals5(1:32,:)<.01);
figure;
X = memTimes(2,:)*1000;
plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to start of pursuit [ms]');ylabel('Number of significant channels (p<.01) /32');
title('MIP Pursuits');legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');title('MIP Pursuit');
%draw vertical line at time=0
line([mean(TT.T66on(trials)-TT.T66on(trials)) mean(TT.T66on(trials)-TT.T66on(trials))],[0 32]);
%draw vertical line at (average) time of reach
line([mean(TT.Treach(trials)-TT.T66on(trials)) mean(TT.Treach(trials)-TT.T66on(trials))],[0 32]);
%draw vertical line when (avg) fix point off, signaling reach
line([mean(TT.T8(trials)-TT.T66on(trials)) mean(TT.T8(trials)-TT.T66on(trials))],[0 32]);
%draw vertical line when direction cue turns on
line([mean(TT.T4(trials)-TT.T66on(trials)) mean(TT.T4(trials)-TT.T66on(trials))],[0 32]);
%draw vertical line when direction cue turns off
line([mean(TT.T6(trials)-TT.T66on(trials)) mean(TT.T6(trials)-TT.T66on(trials))],[0 32]);
axis([-inf inf 0 32]);
%plot figure for PMd channels
sigPVals1 = sum(tuningPVals1(33:48,:)<.01);sigPVals2 = sum(tuningPVals2(33:48,:)<.01);
sigPVals3 = sum(tuningPVals3(33:48,:)<.01);sigPVals4 = sum(tuningPVals4(33:48,:)<.01);
sigPVals5 = sum(tuningPVals5(33:48,:)<.01);
figure;
X = memTimes(2,:)*1000;
plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to start of pursuit [ms]');ylabel('Number of significant channels (p<.01) /16');
title('PMd Pursuits');legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');title('PMd Pursuit');
%draw vertical line at time=0
line([mean(TT.T66on(trials)-TT.T66on(trials)) mean(TT.T66on(trials)-TT.T66on(trials))],[0 16]);
%draw vertical line at (average) time of reach
line([mean(TT.Treach(trials)-TT.T66on(trials)) mean(TT.Treach(trials)-TT.T66on(trials))],[0 16]);
%draw vertical line when (avg) fix point off, signaling reach
line([mean(TT.T8(trials)-TT.T66on(trials)) mean(TT.T8(trials)-TT.T66on(trials))],[0 16]);
%draw vertical line when direction cue turns on
line([mean(TT.T4(trials)-TT.T66on(trials)) mean(TT.T4(trials)-TT.T66on(trials))],[0 16]);
%draw vertical line when direction cue turns off
line([mean(TT.T6(trials)-TT.T66on(trials)) mean(TT.T6(trials)-TT.T66on(trials))],[0 16]);
axis([-inf inf 0 16]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FOR FIXED EYE TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot figure for MIP channels
% sigPVals1 = sum(tuningPVals1(1:32,:)<.01);sigPVals2 = sum(tuningPVals2(1:32,:)<.01);
% sigPVals3 = sum(tuningPVals3(1:32,:)<.01);sigPVals4 = sum(tuningPVals4(1:32,:)<.01);
% sigPVals5 = sum(tuningPVals5(1:32,:)<.01);
% figure;
% X = memTimes(2,:)*1000;
% plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to start of mem period [ms]');ylabel('Number of significant channels (p<.01) /32');
% legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');title('MIP Fixed Eye');
% %draw vertical line at time=0
% line([mean(TT.T6(trials)-TT.T6(trials)) mean(TT.T6(trials)-TT.T6(trials))],[0 32]);
% %draw vertical line at (average) time of reach
% line([mean(TT.Treach(trials)-TT.T6(trials)) mean(TT.Treach(trials)-TT.T6(trials))],[0 32]);
% %draw vertical line when (avg) fix point off, signaling reach
% line([mean(TT.T8(trials)-TT.T6(trials)) mean(TT.T8(trials)-TT.T6(trials))],[0 32]);
% %draw vertical line when direction cue turns on
% line([mean(TT.T4(trials)-TT.T6(trials)) mean(TT.T4(trials)-TT.T6(trials))],[0 32]);
% axis([-inf inf 0 32]);
% %plot figure for PMd channels
% sigPVals1 = sum(tuningPVals1(33:48,:)<.01);sigPVals2 = sum(tuningPVals2(33:48,:)<.01);
% sigPVals3 = sum(tuningPVals3(33:48,:)<.01);sigPVals4 = sum(tuningPVals4(33:48,:)<.01);
% sigPVals5 = sum(tuningPVals5(33:48,:)<.01);
% figure;
% X = memTimes(2,:)*1000;
% plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to start of mem period [ms]');ylabel('Number of significant channels (p<.01) /16');
% legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');title('PMd Fixed Eye');
% %draw vertical line at time=0
% line([mean(TT.T6(trials)-TT.T6(trials)) mean(TT.T6(trials)-TT.T6(trials))],[0 16]);
% %draw vertical line at (average) time of reach
% line([mean(TT.Treach(trials)-TT.T6(trials)) mean(TT.Treach(trials)-TT.T6(trials))],[0 16]);
% %draw vertical line when (avg) fix point off, signaling reach
% line([mean(TT.T8(trials)-TT.T6(trials)) mean(TT.T8(trials)-TT.T6(trials))],[0 16]);
% %draw vertical line when direction cue turns on
% line([mean(TT.T4(trials)-TT.T6(trials)) mean(TT.T4(trials)-TT.T6(trials))],[0 16]);
% axis([-inf inf 0 16]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%