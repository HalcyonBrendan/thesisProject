

%create three time-bin sums in ~15Hz bands

zeroBand = zeros(1,576);
tenBand = zeros(1,576);
twentyBand = zeros(1,576);
thirtyFiveBand = zeros(1,576);
fiftyBand = zeros(1,576);
seventyFiveBand = zeros(1,576);
hundredBand = zeros(1,576);
hundredFiftyBand = zeros(1,576);
twoHundredBand = zeros(1,576);
threeHundredBand = zeros(1,576);

for time = 1:576

    zeroBand(time) = sum(sum(currMag(1:5,time:time+2)));
    tenBand(time) = sum(sum(currMag(6:10,time:time+2)));
    twentyBand(time) = sum(sum(currMag(11:17,time:time+2)));
    thirtyFiveBand(time) = sum(sum(currMag(18:25,time:time+2)));
    fiftyBand(time) = sum(sum(currMag(26:38,time:time+2)));
    seventyFiveBand(time) = sum(sum(currMag(39:51,time:time+2)));
    hundredBand(time) = sum(sum(currMag(52:77,time:time+2)));
    hundredFiftyBand(time) = sum(sum(currMag(78:102,time:time+2)));
    twoHundredBand(time) = sum(sum(currMag(103:154,time:time+2)));
    threeHundredBand(time) = sum(sum(currMag(155:size(chanArray{i,1}.spects,1),time:time+2)));
    
end



%produce spectra for each trial (time independent)

nfft = 2^nextpow2(9501);
nw=4; %default
fs = 1000; %ie. millisecond resolution
for trial = 1:18
   figure;
   pmtm(double(LFP.AD(trial,:)),nw,nfft,fs); 
    
end

%display spectrogram spect
figure;
logSpect = log(squeeze(chanArray{1,1}.spectsAlign(:,:,1)));
imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{1,1}.f(1,:),logSpect);
set(gca,'Ydir','Normal');

for i = 1:8
    figure;
    logSpect = log(mean(chanArray{i,1}.spects,3));
    imagesc(chanArray{i,1}.t(1,:),chanArray{i,1}.f(1,:),logSpect);
    set(gca,'Ydir','Normal');
end

for i = 1:8
    figure;
    logSpect = log(chanArray{i,1}.avgSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),logSpect);
    set(gca,'Ydir','Normal');
end
%look types
for i = 1:8
    figure;
    sacLogSpect = log(chanArray{i,1}.avgSacSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),sacLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/sacAvgCh%d',i);
    savefig(filename);
    figure;
    purLogSpect = log(chanArray{i,1}.avgPurSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),purLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/purAvgCh%d',i);
    savefig(filename);
    figure;
    fixLogSpect = log(chanArray{i,1}.avgFixSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),fixLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/fixAvgCh%d',i);
    savefig(filename);
end
%look dirs
for i = 1:8
    figure;
    zerLogSpect = log(chanArray{i,1}.avgZerSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),zerLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/zerAvgCh%d',i);
    savefig(filename);
    figure;
    ninLogSpect = log(chanArray{i,1}.avgNinSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),ninLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/ninAvgCh%d',i);
    savefig(filename);
    figure;
    oneLogSpect = log(chanArray{i,1}.avgOneSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),oneLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/oneAvgCh%d',i);
    savefig(filename);
    figure;
    twoLogSpect = log(chanArray{i,1}.avgTwoSpect);
    imagesc(1:size(chanArray{i,1}.spects,2)+200,chanArray{i,1}.f(1,:),twoLogSpect);
    set(gca,'Ydir','Normal');
    filename = sprintf('../analysis/mt_spectrograms/twoAvgCh%d',i);
    savefig(filename);
end

%average over different eye movements
sacAvg = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2));
purAvg = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2));
fixAvg = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2));
sacCount = 0;
purCount = 0;
fixCount = 0;

for i = 1:size(chanArray{i,1}.spects,3)

    if strcmp(trialLabels{i},'saccade')
        sacAvg = sacAvg + squeeze(spects(:,:,i));
        sacCount = sacCount+1;
    elseif strcmp(trialLabels{i}, 'pursuit')
        purAvg = purAvg + squeeze(spects(:,:,i));
        purCount = purCount+1;
    else %fixed
        fixAvg = fixAvg + squeeze(spects(:,:,i));
        fixCount = fixCount+1;
    end

end

sacAvg = sacAvg/sacCount;
purAvg = purAvg/purCount;
fixAvg = fixAvg/fixCount;

%average over entire channel
spectAvg = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2));

for i = 1:size(chanArray{i,1}.spects,3)
    spectAvg = spectAvg + squeeze(chanArray{1,1}.spects(:,:,i));
end
spectAvg = spectAvg/size(chanArray{i,1}.spects,3);

%create a list of trial labels
%65=saccade, 66=smooth pursuit, 67=no eye movement
trialLabels = cell(size(chanArray{i,1}.spects,3),2);
for trial = 1:size(chanArray{i,1}.spects,3)
   
     if Events.strobes(trial, 11) == 65
         trialLabels(trial,1) = cellstr('saccade');
     elseif Events.strobes(trial,11) == 66
         trialLabels(trial,1) = cellstr('pursuit');
     else
         trialLabels(trial,1) = cellstr('fixed');
     end

    trialLabels(trial,2) = num2cell(Events.reachAmpAng(trial,2));
end

tuningPVals1 = zeros(48,36);
tuningPVals2 = zeros(48,36);
tuningPVals3 = zeros(48,36);
tuningPVals4 = zeros(48,36);
tuningPVals5 = zeros(48,36);
chanArray = cell(48,1);
for i = 1:48
    
%align trials

%%%%%%%%%% ALIGNMENT PROCEDURE %%%%%%%%%%%%%
%
% 1. For each trial, find time in TT struct that action of interest 
%    starts (ie. pursuit, saccade, etc)
% 2. Shift the times of that trial (in chanArray{i,1}.t) by the above
%    amount so that tAlign = 0 corresponds to start of action
% 3. Round that time to the nearest time found in "calibrator" vector
% 4. Shift spectrogram for that trial to corresponding matrix element
%    chanArray{i,1}.spectsAlign (note that this will make t=0 correspond
%    to cell 221 for current calibrator vector - see below)
%  
%    calibrator = -5.5:.025:5.5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %load unaligned spectrograms and save into chanArray
    handle = sprintf('../analysis/mt_spectrograms/M20100408_256/mt_spects_ch%d.mat',i);
    chanArray{i,1} = load(handle);
    
    %create trialLabels array on first channel
    if i==1
        trialLabels = cell(size(chanArray{i,1}.spects,3),2);
        for trial = 1:size(chanArray{i,1}.spects,3)
             if Events.strobes(trial, 11) == 65
                 trialLabels(trial,1) = cellstr('saccade');
             elseif Events.strobes(trial,11) == 66
                 trialLabels(trial,1) = cellstr('pursuit');
             else
                 trialLabels(trial,1) = cellstr('fixed');
             end
            trialLabels(trial,2) = num2cell(Events.reachAmpAng(trial,2));
        end
    end
    
    %figure out times to align to for each look type
    alignTimes = zeros(size(chanArray{i,1}.spects,3),1);
    for j = 1:size(chanArray{i,1}.spects,3)
        if strcmp(trialLabels{j},'saccade')
            %to saccade
            alignTimes(j,1) = TT.T65sac(j,1);
        elseif strcmp(trialLabels{j}, 'pursuit')
            %to pursuit start
            alignTimes(j,1) = TT.T66on(j,1);
        else %fixed
            %to start of mem period
            alignTimes(j,1) = TT.T6(j,1);
        end
    end


%now shift t matrix by alignTimes vector
%for i = 9:14
    tMS = chanArray{i,1}.t*1000;
    tAligned = tMS-repmat(alignTimes,1,size(chanArray{i,1}.spects,2));
    tAligned = tAligned/1000;
    chanArray{i,1}.tAlign = tAligned;
%end

%now shift spectrograms accordingly
%specStart = zeros(size(chanArray{i,1}.spects,3),1);

%for i = 1:8%48
    chanArray{i,1}.specStart = zeros(size(chanArray{i,1}.spects,3),1);
    for j = 1:size(chanArray{i,1}.spects,3)
        for idx = 1:441
            if abs(chanArray{i,1}.tAlign(j,1)-calibrator(idx)) <= .025/2+.00001
                chanArray{i,1}.specStart(j) = idx;
            end
        end
    end
%end


%spectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,size(chanArray{i,1}.spects,3));
%for i = 1:8%48
    chanArray{i,1}.spectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,size(chanArray{i,1}.spects,3));
    %for j = 1:size(chanArray{i,1}.spects,3)
    %    for k = 1:size(chanArray{i,1}.spects,1)
    %        for l = 1:size(chanArray{i,1}.spects,2)+200
    %            if l < chanArray{i,1}.specStart(j) || l > chanArray{i,1}.specStart(j)+368
    %                chanArray{i,1}.spectsAlign(k,l,j) = 0;
    %            else
    %                chanArray{i,1}.spectsAlign(k,l,j) = chanArray{i,1}.spects(k,l-chanArray{i,1}.specStart(j)+1,j);
    %            end
    %        end
    %    end
    %end
    
    for j = 1:size(chanArray{i,1}.spects,3)
        chanArray{i,1}.spectsAlign(:,chanArray{i,1}.specStart(j):chanArray{i,1}.specStart(j)+size(chanArray{i,1}.spects(:,:,j),2)-1,j) = chanArray{i,1}.spects(:,:,j);
    end
    if mod(i,4)==0
        handle = sprintf('../analysis/workspaces/M20100408_256/chanArray%d-%dAug0315.mat',i-3,i);
        save(handle,'chanArray','-v7.3');
        chanArray = cell(48,1);
    end
    
   disp(i)
end
    
    %run anova for pursuit trials for each channel at various times in mem period
    %to determine directional tuning

    %purs = group==66;
    %fixs = group==67;
    %sacs = group==65;
    %saccades which are farther than median from reach
    %sacThresh = median((TT.T8(sacs)-TT.T65sac(sacs)));
    %isoSacs = (TT.T8-TT.T65sac) > sacThresh;

%     for i = 1:48
%         if mod(i,4)==1
%             handle = sprintf('../analysis/workspaces/chanArray%d-%dAug0315.mat',i,i+3);
%         end
%         clear chanArray;
%         load(handle);
%     for j = 1:length(memTimes)
%         %consider mean power in ~0-10Hz freq band
%         temp = squeeze(mean(chanArray{i,1}.spectsAlign(1:5,221-6+memTimes(1,j),:)));
%         tuningPVals1(i,j) = anovan(temp(fixs),direction(fixs,1),'display', 'off');
%         %consider mean power in ~10-20Hz freq band
%         temp = squeeze(mean(chanArray{i,1}.spectsAlign(6:11,221-6+memTimes(1,j),:)));
%         tuningPVals2(i,j) = anovan(temp(fixs),direction(fixs,1),'display', 'off');
%         %consider mean power in ~20-30Hz freq band
%         temp = squeeze(mean(chanArray{i,1}.spectsAlign(12:16,221-6+memTimes(1,j),:)));
%         tuningPVals3(i,j) = anovan(temp(fixs),direction(fixs,1),'display', 'off');
%         %consider mean power in ~30-40Hz freq band
%         temp = squeeze(mean(chanArray{i,1}.spectsAlign(17:21,221-6+memTimes(1,j),:)));
%         tuningPVals4(i,j) = anovan(temp(fixs),direction(fixs,1),'display', 'off');
%         %consider mean power in ~40-60Hz freq band
%         temp = squeeze(mean(chanArray{i,1}.spectsAlign(22:30,221-6+memTimes(1,j),:)));
%         tuningPVals5(i,j) = anovan(temp(fixs),direction(fixs,1),'display', 'off');
%     end
%     disp(i);
%     end
%     sigPVals1 = sum(tuningPVals1<.05);sigPVals2 = sum(tuningPVals2<.05);
%     sigPVals3 = sum(tuningPVals3<.05);sigPVals4 = sum(tuningPVals4<.05);
%     sigPVals5 = sum(tuningPVals5<.05);
%     figure;
%     plot(X,sigPVals1,X,sigPVals2,X,sigPVals3,X,sigPVals4,X,sigPVals5);xlabel('Time rel. to saccade');ylabel('Number of significant channels (p<.05)');
%     legend('0-10Hz','10-20Hz','20-30Hz','30-40Hz','40-60Hz','Location','northwest');
%     line([mean(TT.T8(fixs)-TT.T6(fixs)) mean(TT.T8(fixs)-TT.T6(fixs))],[0 48]);
%     axis([-inf inf 0 48]);
%     if mod(i,4)==0
%         handle = sprintf('../analysis/workspaces/M20100302_245/chanArray%d-%dAug0315.mat',i-3,i);
%         save(handle,'chanArray','-v7.3');
%         chanArray = cell(48,1);
%     end
%     
%    disp(i)
% end

%take average of spectrograms for look type
for i = 1:8
    chanArray{i,1}.avgSacSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    chanArray{i,1}.avgPurSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    chanArray{i,1}.avgFixSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    for k = 1:size(chanArray{i,1}.spects,1)
        for l = 1:size(chanArray{i,1}.spects,2)+200
            vect1 = chanArray{i,1}.sacSpectsAlign(k,l,:);
            chanArray{i,1}.avgSacSpect(k,l) = mean(vect1(vect1~=0));
            vect2 = chanArray{i,1}.purSpectsAlign(k,l,:);
            chanArray{i,1}.avgPurSpect(k,l) = mean(vect2(vect2~=0));
            vect3 = chanArray{i,1}.fixSpectsAlign(k,l,:);
            chanArray{i,1}.avgFixSpect(k,l) = mean(vect3(vect3~=0));
        end
    end
    disp(i);
end

%take average of spectrograms for look direction
for i = 1:8
    chanArray{i,1}.avgZerSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    chanArray{i,1}.avgNinSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    chanArray{i,1}.avgOneSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    chanArray{i,1}.avgTwoSpect = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200);
    for k = 1:size(chanArray{i,1}.spects,1)
        for l = 1:size(chanArray{i,1}.spects,2)+200
            vect1 = chanArray{i,1}.zerSpectsAlign(k,l,:);
            chanArray{i,1}.avgZerSpect(k,l) = mean(vect1(vect1~=0));
            vect2 = chanArray{i,1}.ninSpectsAlign(k,l,:);
            chanArray{i,1}.avgNinSpect(k,l) = mean(vect2(vect2~=0));
            vect3 = chanArray{i,1}.oneSpectsAlign(k,l,:);
            chanArray{i,1}.avgOneSpect(k,l) = mean(vect3(vect3~=0));
            vect4 = chanArray{i,1}.twoSpectsAlign(k,l,:);
            chanArray{i,1}.avgTwoSpect(k,l) = mean(vect4(vect4~=0));
        end
    end
    disp(i);
end

%get already aligned spectrograms from spectsAlign
for i = 1:8
    sacCount = 1;
    purCount = 1;
    fixCount = 1;
    chanArray{i,1}.sacSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,169);
    chanArray{i,1}.purSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,169);
    chanArray{i,1}.fixSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,108);
    
    for j = 1:size(chanArray{i,1}.spects,3)
        if strcmp(trialLabels{j},'saccade')
            chanArray{i,1}.sacSpectsAlign(:,:,sacCount) = chanArray{i,1}.spectsAlign(:,:,j);
            sacCount = sacCount+1;
        elseif strcmp(trialLabels{j}, 'pursuit')
            chanArray{i,1}.purSpectsAlign(:,:,purCount) = chanArray{i,1}.spectsAlign(:,:,j);
            purCount = purCount+1;
        else
            chanArray{i,1}.fixSpectsAlign(:,:,fixCount) = chanArray{i,1}.spectsAlign(:,:,j);
            fixCount = fixCount+1;
        end
    end
    disp(i);
end

%get already aligned spectrograms from spectsAlign
for i = 1:8
    zerCount = 1;
    ninCount = 1;
    oneCount = 1;
    twoCount = 1;
    chanArray{i,1}.zerSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,113);
    chanArray{i,1}.ninSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,114);
    chanArray{i,1}.oneSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,105);
    chanArray{i,1}.twoSpectsAlign = zeros(size(chanArray{i,1}.spects,1),size(chanArray{i,1}.spects,2)+200,114);
    for j = 1:size(chanArray{i,1}.spects,3)
        if trialLabels{j,2}==0
            chanArray{i,1}.zerSpectsAlign(:,:,zerCount) = chanArray{i,1}.spectsAlign(:,:,j);
            zerCount = zerCount+1;
        elseif trialLabels{j,2}==90
            chanArray{i,1}.ninSpectsAlign(:,:,ninCount) = chanArray{i,1}.spectsAlign(:,:,j);
            ninCount = ninCount+1;
        elseif trialLabels{j,2}==180
            chanArray{i,1}.oneSpectsAlign(:,:,oneCount) = chanArray{i,1}.spectsAlign(:,:,j);
            oneCount = oneCount+1;
        else%==270
            chanArray{i,1}.twoSpectsAlign(:,:,twoCount) = chanArray{i,1}.spectsAlign(:,:,j);
            twoCount = twoCount+1;
        end
    end
    disp(i);
end



