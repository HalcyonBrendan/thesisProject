% %%
% Determine which frequency bands are directionally tuned
% 
% lookType = 'purs'; chanType = 'pmd';
% folders = {'M20100301_246';'M20100407_456';'M20100408_256'};
% %folders = {'M20100302_645';'M20100302_245'};
% sigBandCounts = zeros(13,2);
% 
% for j = 1:2
%     if j == 1
%         trainType = 'cis';
%     else
%         trainType = 'trans';
%     end
%     numChans = 0;
%     for i = 1:length(folders)
%         folder = folders{i};
%         % Get tuning information
%         handle = sprintf('../analysis/workspaces/%s/%s_dirANOVA_240416.mat',folder,trainType);
%         tunData = load(handle);
%         % Set p-value you want to look at
%         p = 0.05;
% 
%         % Get the appropriate info
%         if strcmp(lookType,'sacs')
%             tuningPVals = tunData.sacTuningPVals;
%             if strcmp(trainType,'cis')
%                 trials = tunData.cisSacs;
%             elseif strcmp(trainType,'trans')
%                 trials = tunData.transSacs;
%             end
%         elseif strcmp(lookType,'purs')
%             tuningPVals = tunData.purTuningPVals;
%             if strcmp(trainType,'cis')
%                 trials = tunData.cisPurs;
%             elseif strcmp(trainType,'trans')
%                 trials = tunData.transPurs;
%             end
%         end
% 
% 
%         if strcmp(chanType,'mip')
%             numChans = numChans + 32;
%             sigBands = s(tuningPVals(1:32,1,:,4)<p);
%         elseif strcmp(chanType,'pmd')
%             numChans = numChans + 16;
%             sigBands = s(tuningPVals(33:48,1,:,4)<p);
%         end
% 
%         sigBandCounts(:,j) = sigBandCounts(:,j) + sum(sigBands,1)';
%     end
% end
% sigBandCounts = sigBandCounts/numChans;
% 
% % Produce figure
% figure; b=bar(sigBandCounts,'BarWidth',1.3);
% set(b(1),'FaceColor','r');%set(b(2),'FaceColor','r');
% ax = gca; ax.XTickLabels = ['  0-5 Hz  ';' 5-15 Hz  ';' 15-25 Hz ';' 25-35 Hz ';' 35-45 Hz ';' 45-55 Hz ';' 55-65 Hz ';' 65-75 Hz ';' 75-85 Hz ';' 85-95 Hz ';' 95-105 Hz';'105-125 Hz';'125-150 Hz'];
% ax.XTickLabelRotation=60; ylim([0 .63]);
% fig=gcf; fig.Position = [150 150 300 200];
% 
% 
% %%
% % Get direction classification results for single data sets
% trainType = 'cis'; lookType = 'sacs'; chanType = 'mip'; numChans = 12;
% %folders = {'M20100301_246';'M20100407_456';'M20100408_256','M20100302_645';'M20100302_245'};
% folders = {'M20100301_246'};
% numFolders = length(folders);
% 
% for fold = 1:numFolders
%     
% folder = folders{fold};
% handle = sprintf('../data/%s/LFPChan01.mat',folder);
% load(handle);
% handle = sprintf('../analysis/workspaces/%s/%s%s_dirClass_%dchs_240416_%s.mat',folder,lookType,chanType,numChans,trainType);
% load(handle);
% info = Events.dataSummary;
% 
% numIts = size(valTestTrials,2);
% 
% numCis = 0; numTrans = 0; numTs = size(meanPercentCorr,1);
% cisCorr = zeros(1,numTs); transCorr = zeros(1,numTs);
% cisAns = zeros(100000,numTs); transAns = zeros(100000,numTs);
% cIdx = 1; tIdx = 1;
% % Get results for each type of look
% for i = 1:numIts
%     % Find cis and trans indices
%     cisTests = find(ismember(valTestTrials(:,i),valCisTrials));
%     transTests = find(ismember(valTestTrials(:,i),valTransTrials));
%     
%     numCurrCis = length(cisTests); numCurrTrans = length(transTests);
%     numCis = numCis + numCurrCis;
%     numTrans = numTrans + numCurrTrans;
%     
%     cisAns(cIdx:cIdx+numCurrCis-1,:) = s(answerCheck(:,i,1,cisTests,3))';
%     transAns(tIdx:tIdx+numCurrTrans-1,:) = s(answerCheck(:,i,1,transTests,3))';
%     cIdx = cIdx + numCurrCis; tIdx = tIdx + numCurrTrans;
%     
%     cisCorr = cisCorr + s(sum(answerCheck(:,i,1,cisTests,3),4))';
%     transCorr = transCorr + s(sum(answerCheck(:,i,1,transTests,3),4))';
% 
% end
% 
% cisPer = cisCorr/numCis;
% transPer = transCorr/numTrans;
% 
% % Compute bootstrap confidence intervals
% numBoots = 100000;
% if max(cisAns(:,1))>0
%     numTests = length(valCisTrials);
%     cisBootNumCorr = zeros(numBoots,numTs);
%     for i = 1:numBoots
%         % Select indices (with replacement) to sample
%         bootIdxs = randperm(numCis,numTests);
%         cisBootNumCorr(i,:) = sum(cisAns(bootIdxs,:),1);
%     end
%     cisBootNumCorr = sort(cisBootNumCorr)/numTests;
%     cisCI95LB = cisBootNumCorr(round(.05*numBoots),:);
%     cisCI95UB = cisBootNumCorr(round(.95*numBoots),:);
%     cisCI99LB = cisBootNumCorr(round(.01*numBoots),:);
%     cisCI99UB = cisBootNumCorr(round(.99*numBoots),:);
% end
% if sum(transAns(:,1))>0
%     numTests = length(valTransTrials)+1;
%     transBootNumCorr = zeros(numBoots,numTs);
%     for i = 1:numBoots
%         % Select indices (with replacement) to sample
%         bootIdxs = randperm(numTrans,numTests);
%         transBootNumCorr(i,:) = sum(transAns(bootIdxs,:),1);
%     end
%     transBootNumCorr = sort(transBootNumCorr)/numTests;
%     transCI95LB = transBootNumCorr(round(.05*numBoots),:);
%     transCI95UB = transBootNumCorr(round(.95*numBoots),:);
%     transCI99LB = transBootNumCorr(round(.01*numBoots),:);
%     transCI99UB = transBootNumCorr(round(.99*numBoots),:);
% end
% 
% handle = sprintf('../analysis/workspaces/%s/%s%s_dirResults_%dchs_240416_%s.mat',folder,type,chanType,numChans,dirType);
% if strcmp(folder,'M20100302_645') || strcmp(folder,'M20100302_245')
%     save(handle,'cisBootNumCorr','cisCI95LB','cisCI95UB','cisCI99LB','cisCI99UB','numBoots','numTests','cisCorr','transCorr','cisPer','transPer','numCis','numTrans','trainType','lookType','chanType','answerCheck','channels','folder','meanPercentCorr','numCorr','percentCorr','timeIdxs','valTestTrials','valTrainTrials','valCisTrials','valTransTrials');
% else
%     save(handle,'cisBootNumCorr','transBootNumCorr','cisCI95LB','cisCI95UB','cisCI99LB','cisCI99UB','transCI95LB','transCI95UB','transCI99LB','transCI99UB','numBoots','numTests','cisCorr','transCorr','cisPer','transPer','numCis','numTrans','trainType','lookType','chanType','answerCheck','channels','folder','meanPercentCorr','numCorr','percentCorr','timeIdxs','valTestTrials','valTrainTrials','valCisTrials','valTransTrials');
% end
% %clearvars -except numChans trainType lookType chanType folders numFolders fold
% 
% end
% 
% %%
% % Get direction classification results for pooled data cis/trans data sets
% 
% chType = 'pmd';
% lookType = 'purs';
% trainType = 'trans';
% 
% folders = {'M20100301_246';'M20100407_456'};%;'M20100408_256'};
% 
% numCis = 0; numTrans = 0; numCisTests = 0; numTransTests = 0; numTs = 4; numChans = 12;
% cisCorr = zeros(1,numTs); transCorr = zeros(1,4); cIdx = 1; tIdx = 1;
% cisAns = zeros(100000,numTs); transAns = zeros(100000,numTs);
% for j = 1:length(folders)
%     handle = sprintf('../analysis/workspaces/%s/%s%s_dirClass_12chs_240416_%s.mat',folders{j},lookType,chType,trainType);
%     load(handle);
%     handle = sprintf('../data/%s/LFPChan01.mat',folders{j});
%     load(handle);
%     info = Events.dataSummary;
%     
%     numIts = size(valTestTrials,2);
%     numCisTests = numCisTests+length(valCisTrials);
%     numTransTests = numTransTests+length(valTransTrials);
%     
%     for i = 1:numIts        
%         % Find cis and trans indices
%         cisTests = find(ismember(valTestTrials(:,i),valCisTrials));
%         transTests = find(ismember(valTestTrials(:,i),valTransTrials));
% 
%         numCurrCis = length(cisTests); numCurrTrans = length(transTests);
%         numCis = numCis + numCurrCis; numTrans = numTrans + numCurrTrans;
% 
%         cisAns(cIdx:cIdx+numCurrCis-1,:) = s(answerCheck(:,i,1,cisTests,3))';
%         transAns(tIdx:tIdx+numCurrTrans-1,:) = s(answerCheck(:,i,1,transTests,3))';
%         cIdx = cIdx + numCurrCis; tIdx = tIdx + numCurrTrans;
% 
%         cisCorr = cisCorr + s(sum(answerCheck(:,i,1,cisTests,3),4))';
%         transCorr = transCorr + s(sum(answerCheck(:,i,1,transTests,3),4))';
%     end
% end
% 
% cisPer = cisCorr/numCis;
% transPer = transCorr/numTrans;
% 
% % Compute bootstrap confidence intervals
% numBoots = 100000;
% if max(cisAns(:,1))>0
%     cisBootNumCorr = zeros(numBoots,numTs);
%     for i = 1:numBoots
%         % Select indices (with replacement) to sample
%         bootIdxs = randperm(numCis,numCisTests);
%         cisBootNumCorr(i,:) = sum(cisAns(bootIdxs,:),1);
%     end
%     cisBootNumCorr = sort(cisBootNumCorr)/numCisTests;
%     cisCI95LB = cisBootNumCorr(round(.05*numBoots),:);
%     cisCI95UB = cisBootNumCorr(round(.95*numBoots),:);
%     cisCI99LB = cisBootNumCorr(round(.01*numBoots),:);
%     cisCI99UB = cisBootNumCorr(round(.99*numBoots),:);
% end
% if max(transAns(:,1))>0
%     transBootNumCorr = zeros(numBoots,numTs);
%     for i = 1:numBoots
%         % Select indices (with replacement) to sample
%         bootIdxs = randperm(numTrans,numTransTests);
%         transBootNumCorr(i,:) = sum(transAns(bootIdxs,:),1);
%     end
%     transBootNumCorr = sort(transBootNumCorr)/numTransTests;
%     transCI95LB = transBootNumCorr(round(.05*numBoots),:);
%     transCI95UB = transBootNumCorr(round(.95*numBoots),:);
%     transCI99LB = transBootNumCorr(round(.01*numBoots),:);
%     transCI99UB = transBootNumCorr(round(.99*numBoots),:);
% end
% 
% handle = sprintf('../analysis/workspaces/%s%s_pooledDirResults_%dchs_240416_%sBoth.mat',type,chanType,numChans,dirType);
% save(handle,'cisBootNumCorr','cisCI95LB','cisCI95UB','cisCI99LB','cisCI99UB','numBoots','numCisTests','numTransTests','cisCorr','cisPer','numCis','trainType','lookType','chanType','channels','timeIdxs','transBootNumCorr','transCI95LB','transCI95UB','transCI99LB','transCI99UB','transCorr','transPer','numTrans');

%%%
% % Get direction classification results for cis only data sets
% 
% chType = 'pmd';
% lookType = 'purs';
% trainType = 'cis';
% 
% folders = {'M20100302_645';'M20100302_245'};
% 
% numCis = 0; numTests = 0; numTs = 4; numChans = 12;
% cisCorr = zeros(1,numTs); transCorr = zeros(1,4); cIdx = 1;
% cisAns = zeros(100000,numTs);
% for j = 1:length(folders)
%     handle = sprintf('../analysis/workspaces/%s/%s%s_dirClass_12chs_240416_%s.mat',folders{j},lookType,chType,trainType);
%     load(handle);
%     handle = sprintf('../data/%s/LFPChan01.mat',folders{j});
%     load(handle);
%     info = Events.dataSummary;
%     
%     numIts = size(valTestTrials,2);
%     numTests = numTests+length(valCisTrials);
%     
%     for i = 1:numIts        
%         % Find cis indices
%         cisTests = find(ismember(valTestTrials(:,i),valCisTrials));
% 
%         numCurrCis = length(cisTests);
%         numCis = numCis + numCurrCis;
% 
%         cisAns(cIdx:cIdx+numCurrCis-1,:) = s(answerCheck(:,i,1,cisTests,3))';
%         cIdx = cIdx + numCurrCis;
% 
%         cisCorr = cisCorr + s(sum(answerCheck(:,i,1,cisTests,3),4))';
%     end
% end
% 
% cisPer = cisCorr/numCis;
% 
% % Compute bootstrap confidence intervals
% numBoots = 100000;
% if max(cisAns(:,1))>0
%     cisBootNumCorr = zeros(numBoots,numTs);
%     for i = 1:numBoots
%         % Select indices (with replacement) to sample
%         bootIdxs = randperm(numCis,numTests);
%         cisBootNumCorr(i,:) = sum(cisAns(bootIdxs,:),1);
%     end
%     cisBootNumCorr = sort(cisBootNumCorr)/numTests;
%     cisCI95LB = cisBootNumCorr(round(.05*numBoots),:);
%     cisCI95UB = cisBootNumCorr(round(.95*numBoots),:);
%     cisCI99LB = cisBootNumCorr(round(.01*numBoots),:);
%     cisCI99UB = cisBootNumCorr(round(.99*numBoots),:);
% end
% 
% handle = sprintf('../analysis/workspaces/%s%s_pooledDirResults_%dchs_240416_%sOnly.mat',type,chanType,numChans,dirType);
% save(handle,'cisBootNumCorr','cisCI95LB','cisCI95UB','cisCI99LB','cisCI99UB','numBoots','numTests','cisCorr','cisPer','numCis','trainType','lookType','chanType','channels','timeIdxs');
%
%
% %%
%
% Make plots for classification rates at each time period
% mip=load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/sacsmip_pooledDirResults_12chs_240416_cisOnly.mat');
% pmd=load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/sacspmd_pooledDirResults_12chs_240416_cisOnly.mat');
% 
% % Cis only:
% corrPer = zeros(1,9); lb = zeros(1,9); ub = zeros(1,9);
% corrPer(5)=NaN; lb(5)=NaN; ub(5)=NaN;
% % MIP
% % Compute size of error bars
% lb(1:4) = mip.cisPer-mip.cisCI95LB;
% ub(1:4) = mip.cisCI95UB-mip.cisPer;
% corrPer(1:4) = mip.cisPer;
% 
% % PMd
% % Compute size of error bars
% lb(6:9) = pmd.cisPer-pmd.cisCI95LB;
% ub(6:9) = pmd.cisCI95UB-pmd.cisPer;
% corrPer(6:9) = pmd.cisPer;
% 
% % Add to plot
% for i = 1:9;
%     if i < 5
%         col='g';
%     else
%         col='c';
%     end
%     bar(i,corrPer(i),'BarWidth',0.8,'FaceColor',col); hold on;
%     errorbar(i,corrPer(i),lb(i),ub(i),'ko');
% end
% line([0 10],[.25 .25],'color','k','LineWidth',1.7,'LineStyle','--');
% set(gca, 'XTick', [1:4 6:9]);
% set(gca, 'XTickLabel', {'T1' 'T2' 'T3' 'T4' 'T1' 'T2' 'T3' 'T4'});
% line([5 5],[0 .5],'color','k','LineWidth',1);
% fig =gcf; fig.Position = [150 150 570 210];
% ylim([0 .5]);
% % %%
% % Make plots for classification rates at each time period
% mip=load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/pursmip_pooledDirResults_12chs_240416_transBoth.mat');
% pmd=load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/purspmd_pooledDirResults_12chs_240416_transBoth.mat');
% 
% corr1 = mip.cisPer;
% corr1 = [corr1; mip.transPer]';
% corr2 = pmd.cisPer;
% corr2 = [corr2; pmd.transPer]';
% lb1 = mip.cisPer - mip.cisCI95LB;
% lb1 = [lb1; mip.transPer - mip.transCI95LB];
% ub1 = abs(mip.cisPer - mip.cisCI95UB);
% ub1 = [ub1; abs(mip.transPer - mip.transCI95UB)];
% lb2 = pmd.cisPer - pmd.cisCI95LB;
% lb2 = [lb2; pmd.transPer - pmd.transCI95LB];
% ub2 = abs(pmd.cisPer - pmd.cisCI95UB);
% ub2 = [ub2; abs(pmd.transPer - pmd.transCI95UB)];
% 
% figure;
% bar1=bar(1:3,corr1(2:4,:),'BarWidth',1); hold on;
% set(bar1(1),'FaceColor','g');set(bar1(2),'FaceColor','b');
% bar2=bar(5:7,corr2(2:4,:),'BarWidth',1);
% set(bar2(1),'FaceColor','c');set(bar2(2),'FaceColor','m');
% errorbar(.85:2.85,corr1(2:4,1),lb1(1,2:4),ub1(1,2:4),'ko');
% errorbar(1.15:3.15,corr1(2:4,2),lb1(2,2:4),ub1(2,2:4),'ko');
% errorbar(4.85:6.85,corr2(2:4,1),lb2(1,2:4),ub2(1,2:4),'ko');
% errorbar(5.15:7.15,corr2(2:4,2),lb2(2,2:4),ub2(2,2:4),'ko');
% 
% xlim([0 8]);ylim([0 .55]);
% line([0 10],[.25 .25],'color','k','LineWidth',1.7,'LineStyle','--');
% set(gca, 'XTick', [1:3 5:7]);
% set(gca, 'XTickLabel', {'T2' 'T3' 'T4' 'T2' 'T3' 'T4'});
% line([4 4],[0 .6],'color','k','LineWidth',1);
% fig =gcf; fig.Position = [150 150 570 240];

% %%
% Compute percentage of correct classifications for moves opposite to
% reaches vs perpendicular to reaches

lookType = 'purs'; trainType = 'trans';
folders = {'M20100301_246';'M20100407_456'};%'M20100408_256'};
%folders = {'M20100302_645';'M20100302_245'};
for j = 1:2
    if j == 1
        chanType = 'mip';
    else
        chanType = 'pmd';
    end
    
    numOpps = 0; numPerps = 0; numOppCorr = zeros(4,1); numPerpCorr = zeros(4,1);
    for k = 1:length(folders)
        
        folder = folders{k};
        handle = sprintf('../data/%s/LFPChan01.mat',folder);
        load(handle);
        info = Events.dataSummary;
        handle = sprintf('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/%s/%s%s_dirResults_12chs_240416_%s.mat',folder,lookType,chanType,trainType);
        data = load(handle);
        
        eyeDir = zeros(length(info),1);
        eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==11) = 1;
        eyeDir(Events.EyeHandCoord(:,1)==0 & Events.EyeHandCoord(:,2)==-9) = 2;
        eyeDir(Events.EyeHandCoord(:,1)==-10 & Events.EyeHandCoord(:,2)==1) = 3;
        eyeDir(Events.EyeHandCoord(:,1)==10 & Events.EyeHandCoord(:,2)==1) = 4;
        
        reachDir = zeros(length(info),1);
        reachDir(Events.EyeHandCoord(:,7)==0 & Events.EyeHandCoord(:,8)==7) = 1;
        reachDir(Events.EyeHandCoord(:,7)==0 & Events.EyeHandCoord(:,8)==-7) = 2;
        reachDir(Events.EyeHandCoord(:,7)==-7 & Events.EyeHandCoord(:,8)==0) = 3;
        reachDir(Events.EyeHandCoord(:,7)==7 & Events.EyeHandCoord(:,8)==0) = 4;
        
        
        numIts = size(data.valTestTrials,2);
        for i = 1:numIts
            trans = data.valTestTrials(ismember(data.valTestTrials(:,i),data.valTransTrials),i);
            cis = data.valTestTrials(ismember(data.valTestTrials(:,i),data.valCisTrials),i);
            
            transLook = eyeDir(trans);
            transReach = reachDir(trans);
            
            oppTrans = (transLook<3 & transReach<3) | (transLook>2 & transReach>2);
            perpTrans = oppTrans==0;
            oppTrans = trans(oppTrans);
            perpTrans = trans(perpTrans);
            numOpps = numOpps +length(oppTrans);
            numPerps = numPerps+length(perpTrans);
            
            transSuccess = zeros(length(trans),4);
            for t = 1:4
                transSuccess(:,t) = ismember(trans,data.valTestTrials(s(data.answerCheck(t,i,1,:,3))==1,i));
                %cisSuccess(:,t) = cis(ismember(cis,data.valTestTrials(s(data.answerCheck(t,i,1,:,3))==1,i)));
                numOppCorr(t) = numOppCorr(t) + sum(ismember(trans(find(transSuccess(:,t))),oppTrans));
                numPerpCorr(t) = numPerpCorr(t) + sum(ismember(trans(find(transSuccess(:,t))),perpTrans));
            end
        end
        
    end
    if j == 1
        mipOppCor = numOppCorr;
        mipOpps = numOpps;
        mipPerpCor = numPerpCorr;
        mipPerps = numPerps;
    else
        pmdOppCor = numOppCorr;
        pmdOpps = numOpps;
        pmdPerpCor = numPerpCorr;
        pmdPerps = numPerps;
    end
end
mipOppPer = mipOppCor/mipOpps;
mipPerpPer = mipPerpCor/mipPerps;
pmdOppPer = pmdOppCor/pmdOpps;
pmdPerpPer = pmdPerpCor/pmdPerps;

% Compute bootstrap confidence intervals
numMOs = ceil(mipOpps/numIts); numMPs = ceil(mipPerps/numIts);
numPOs = ceil(pmdOpps/numIts); numPPs = ceil(pmdPerps/numIts);
meanMO = ceil(mipOppPer(2:4)*1000);
meanMP = ceil(mipPerpPer(2:4)*1000);
meanPO = ceil(pmdOppPer(2:4)*1000);
meanPP = ceil(pmdPerpPer(2:4)*1000);
MOTrials = zeros(1000,3);
MPTrials = zeros(1000,3);
POTrials = zeros(1000,3);
PPTrials = zeros(1000,3);
for t = 1:3
    MOTrials(1:meanMO(t),t) = 1;
    MPTrials(1:meanMP(t),t) = 1;
    POTrials(1:meanPO(t),t) = 1;
    PPTrials(1:meanPP(t),t) = 1;
end
numBoots = 100000;
moBoot = zeros(numBoots,3);mpBoot = zeros(numBoots,3);
poBoot = zeros(numBoots,3);ppBoot = zeros(numBoots,3);
for i = 1:numBoots
    for t = 1:3
        moIdxs = randperm(1000,numMOs);
        mpIdxs = randperm(1000,numMPs);
        poIdxs = randperm(1000,numPOs);
        ppIdxs = randperm(1000,numPPs);
        moBoot(i,t) = sum(MOTrials(moIdxs,t),1);
        mpBoot(i,t) = sum(MPTrials(mpIdxs,t),1);
        poBoot(i,t) = sum(POTrials(poIdxs,t),1);
        ppBoot(i,t) = sum(PPTrials(ppIdxs,t),1);
    end
end

moBoot = sort(moBoot,1);
mpBoot = sort(mpBoot,1);
poBoot = sort(poBoot,1);
ppBoot = sort(ppBoot,1);

moLB95 = moBoot(5000,:)/numMOs;
moUB95 = moBoot(95000,:)/numMOs;
mpLB95 = mpBoot(5000,:)/numMPs;
mpUB95 = mpBoot(95000,:)/numMPs;
poLB95 = poBoot(5000,:)/numPOs;
poUB95 = poBoot(95000,:)/numPOs;
ppLB95 = ppBoot(5000,:)/numPPs;
ppUB95 = ppBoot(95000,:)/numPPs;

moLBSize = mipOppPer(2:4)'-moLB95;
moUBSize = moUB95-mipOppPer(2:4)';
mpLBSize = mipPerpPer(2:4)'-mpLB95;
mpUBSize = mpUB95-mipPerpPer(2:4)';
poLBSize = pmdOppPer(2:4)'-poLB95;
poUBSize = poUB95-pmdOppPer(2:4)';
ppLBSize = pmdPerpPer(2:4)'-ppLB95;
ppUBSize = ppUB95-pmdPerpPer(2:4)';

% Produce bar graph
figure;
bar1=bar(1:3,[mipOppPer(2:4,1) mipPerpPer(2:4,1)],'BarWidth',1); hold on;
set(bar1(1),'FaceColor','g');
set(gca, 'XTickLabel', {'T2' 'T3' 'T4'});
line([0.5 3.5],[.25 .25],'color','k','LineWidth',1.7,'LineStyle','--');
errorbar(.85:2.85,mipOppPer(2:4),moLBSize,moUBSize,'ko');
errorbar(1.15:3.15,mipPerpPer(2:4),mpLBSize,mpUBSize,'ko');
ylim([0 .7]);

% Produce bar graph
figure;
bar1=bar(1:3,[pmdOppPer(2:4,1) pmdPerpPer(2:4,1)],'BarWidth',1); hold on;
set(bar1(1),'FaceColor','c');
set(gca, 'XTickLabel', {'T2' 'T3' 'T4'});
line([0.5 3.5],[.25 .25],'color','k','LineWidth',1.7,'LineStyle','--');
errorbar(.85:2.85,pmdOppPer(2:4),poLBSize,poUBSize,'ko');
errorbar(1.15:3.15,pmdPerpPer(2:4),ppLBSize,ppUBSize,'ko');
ylim([0 .7]);