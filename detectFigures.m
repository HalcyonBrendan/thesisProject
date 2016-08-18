% %%
% % Brendan Thorn
% % April 19, 2016
% %
% % detectFigures.m
% %
% %%
% Build bar graph of positive detects (true or false) vs time for single dataset
%
% 
% chType = 'mip';
% lookType = 'sacs';
% 
% folder = 'M20100302_645';
% handle = sprintf('../analysis/workspaces/%s/%s%s_optRun_12chs_240416.mat',folder,lookType,chType);
% load(handle);
% 
% % Add up positives at each time step
% posByTime = sum(adjResults,1);
% % Pool them into larger time bins
% posByBins = reshape(posByTime(1:40),4,10);
% posByBins = sum(posByBins,1); posByBins(10) = posByBins(10)+posByTime(41);
% Z = sum(posByBins);
% posByBins = posByBins/Z;
% posAfter = zeros(length(posByTime),1);
% posAfter(2:end,1) = posByTime(1:end-1);
% 
% %figure;
% %bar(t(testTimes)*1000+150,posByTime);
% figure;
% bar(-850:100:-50,posByBins(1:9),'FaceColor','r','BarWidth',1.0); hold on;
% bar(50,posByBins(10),'FaceColor','g','BarWidth',100); hold on;
% xlim([-900 100]); ylim([0 .8]);
% ax = gca; ax.XTick = -900:100:100;
% line([0 0],[0 .81],'color','k','LineWidth',3);
% line([-900 100],[0 1],'color','k','LineWidth',2,'LineStyle','--');
% plot(linspace(-900,100,41),cumsum(posAfter)/sum(posAfter),'LineWidth',2.0);

% %%
%Build bar graph of positive detects (true or false) vs time for pooled data
% 
% 
% chType = 'pmd';
% lookType = 'purs';
% 
% folders = {'M20100302_645';'M20100301_246';'M20100302_245';'M20100407_456';'M20100408_256'};
% 
% posByTime = zeros(1,41);
% for i = 1:5
%     
%     handle = sprintf('../analysis/workspaces/%s/%s%s_optRun_12chs_240416.mat',folders{i},lookType,chType);
%     load(handle);
%     
%     % Add up positives at each time step
%     posByTime = posByTime + sum(adjResults,1);
% 
% end
% % Pool them into larger time bins
% posByBinsCurr = reshape(posByTime(1:40),4,10);
% posByBins = sum(posByBinsCurr,1); posByBins(10) = posByBins(10)+posByTime(41);
% 
% posByBins = posByBins/sum(posByBins);
% 
% posAfter = zeros(length(posByTime),1);
% posAfter(2:end,1) = posByTime(1:end-1);
% 
% figure;
% bar(-850:100:-50,posByBins(1:9),'FaceColor','r','BarWidth',1.0); hold on;
% bar(50,posByBins(10),'FaceColor','g','BarWidth',100); hold on;
% xlim([-900 100]); ylim([0 .8]);
% ax = gca; ax.XTick = -900:100:100;
% line([0 0],[0 .81],'color','k','LineWidth',3);
% line([-900 100],[0 1],'color','k','LineWidth',2,'LineStyle','--');
% plot(linspace(-900,100,41),cumsum(posAfter)/sum(posAfter),'LineWidth',2.0);

% %%
% % Build detects plot like in svm_lookDetect3*, but with solid exec period
% 
% % Take max of execution period to check for true detection
% % Multiply to get specific colors
% plotResults = .8*adjResults;
% plotResults(:,37:41) = .5*repmat(max(plotResults(:,37:41),[],2),1,5);
% plotResults(plotResults(:,37)==0,37:41) = 0.55;
% plotResults(end,35)=1;
% 
% figure;imagesc(1000*t(testTimes)+155,1:numTestTrials,plotResults);colormap(jet(8));%colorbar();
% ax=gca;ax.YDir='normal';ax.XTickLabel = -900:100:100;
% %xlabel('Time Relative to Saccade [ms]');ylabel('Test Trial');
% %tit = sprintf('%s %s Detects SVM 170416',chType,lookType);title(tit); 
% xlim([-900 100]);
% line([-5 -5],[0 numTestTrials+1],'color','k','LineWidth',4);
% %line([-.15 -.15],[0 numTestTrials],'Color','r','LineWidth',2);
% % %line([0 0],[0 numTestTrials],'Color','g','LineWidth',2);line([-.05 -.05],[0 numTestTrials],'Color','c','LineWidth',2);

% %%
% % Make simple bar graph with error bars of false positive/true positive for all datasets
% 
% chType = 'pmd';
% lookType = 'purs';
% 
% folders = {'M20100302_645';'M20100301_246';'M20100302_245';'M20100407_456';'M20100408_256'};
% numSets = length(folders);
% 
% movesCorr = zeros(numSets,1);
% lB = zeros(numSets,1);
% uB = zeros(numSets,1);
% lBsize = zeros(numSets,1);
% uBsize = zeros(numSets,1);
% totTrials = 0; totNumCorr = 0;
% for i = 1:numSets
%     % Load results
%     handle = sprintf('../analysis/workspaces/%s/%s%s_optRun_12chs_240416.mat',folders{i},lookType,chType);
%     load(handle);
%     
%     % Get percentage of correct moves
%     movesCorr(i) = movesCorrPerc;
%     numCorr = round(movesCorrPerc*numTestTrials);
%     totNumCorr = totNumCorr + numCorr;
%     totTrials = totTrials + numTestTrials;
%     % Get lengths of confidence intervals
%     lB(i) = 1-betainv(.975,numTestTrials-numCorr+1,numCorr);
%     uB(i) = 1-betainv(.025,numTestTrials-numCorr,numCorr+1);
%     
%     lBsize(i) = movesCorrPerc-lB(i);
%     uBsize(i) = uB(i)-movesCorrPerc;
% end
% 
% pooledCorr = totNumCorr/totTrials;
% pooledLB = 1-betainv(.975,totTrials-totNumCorr+1,totNumCorr);
% pooledUB = 1-betainv(.025,totTrials-totNumCorr,totNumCorr+1);
% pooledLBsize = pooledCorr-pooledLB;
% pooledUBsize = pooledUB-pooledCorr;

% % Make plot
% figure;
% bar(1:5,movesCorr,'FaceColor','c','BarWidth',0.8); hold on;
% bar(6.5,pooledCorr,'FaceColor','m','BarWidth',0.8); hold on;
% errorbar(1:5,movesCorr,lBsize,uBsize,'ko'); hold on;
% errorbar(6.5,pooledCorr,pooledLBsize,pooledUBsize,'ko'); hold off;
% line([0 7],[.226 .226],'color','b','LineWidth',1.7);
% line([5.75 5.75],[0 0.7],'color','k','LineWidth',2);
% ax=gca; ax.XTick = [1:5 6.5];
% ax.XTickLabel = [' M645 ';' M246 ';' M245 ';' M456 ';' M256 ';'Pooled'];
% %%
% Make plot of only pooled values
% figure;
% bar(1,corr(1),'FaceColor','g','BarWidth',0.8); hold on;
% bar(2,corr(2),'FaceColor','c','BarWidth',0.8); hold on;
% bar(3,corr(3),'FaceColor','m','BarWidth',0.8); hold on;
% bar(4,corr(4),'FaceColor','r','BarWidth',0.8); hold on;
% errorbar(1,corr(1),lb(1),ub(1),'ko'); hold on;
% errorbar(2,corr(2),lb(2),ub(2),'ko'); hold on;
% errorbar(3,corr(3),lb(3),ub(3),'ko'); hold on;
% errorbar(4,corr(4),lb(4),ub(4),'ko');
% line([0 5],[.226 .226],'color','k','LineWidth',1.5,'LineStyle','--');
% ax=gca; ax.XTick = 1:4;
% ax.XTickLabel = ['MIP Saccades'; 'PMd Saccades'; 'MIP Pursuits'; 'PMd Pursuits'];

% %%
% Determine which frequency bands are temporally tuned

lookType = 'purs'; chanType = 'pmd';
folders = {'M20100301_246';'M20100407_456';'M20100408_256';'M20100302_645';'M20100302_245'};
sigBandCounts = zeros(13,1);

numChans = 0;
for i = 1:length(folders)
    folder = folders{i};
    % Get tuning information
    handle = sprintf('../analysis/workspaces/%s/detANOVA_170416_6973v7781.mat',folder);
    tunData = load(handle);
    % Set p-value you want to look at
    p = 0.001;

    % Get the appropriate info
    if strcmp(lookType,'sacs')
        tuningPVals = tunData.sacDetectPVals;
    elseif strcmp(lookType,'purs')
        tuningPVals = tunData.purDetectPVals;
    end

    if strcmp(chanType,'mip')
        numChans = numChans + 32;
        sigBands = s(tuningPVals(1:32,1,:)<p);
    elseif strcmp(chanType,'pmd')
        numChans = numChans + 16;
        sigBands = s(tuningPVals(33:48,1,:)<p);
    end

    sigBandCounts = sigBandCounts + sum(sigBands,1)';
end

sigBandCounts = sigBandCounts/numChans;

% Produce figure
figure; b=bar(sigBandCounts,'BarWidth',.9);
set(b(1),'FaceColor','r');%set(b(2),'FaceColor','r');
ax = gca; ax.XTickLabels = ['  0-5 Hz  ';' 5-15 Hz  ';' 15-25 Hz ';' 25-35 Hz ';' 35-45 Hz ';' 45-55 Hz ';' 55-65 Hz ';' 65-75 Hz ';' 75-85 Hz ';' 85-95 Hz ';' 95-105 Hz';'105-125 Hz';'125-150 Hz'];
ax.XTickLabelRotation=35;
ylim([0 0.6]);
fig =gcf; fig.Position = [150 150 570 225];

