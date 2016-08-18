[freqBands, fBounds] = create_freqBands(spects, f, fBounds);
% size(freqBands)
% ans = 32    13   140   446

% First choice: ch1, b12, p=2.8e-7
% Second choice: ch17, b5, p=0.040

pows = s(freqBands(17,5,81,:));
info = Events.dataSummary;

figure;
plot(pows(info(:,1)==129 & info(:,4)==65,1),'rx');hold on;
plot(pows(info(:,1)==130 & info(:,4)==65,1),'m^');
plot(pows(info(:,1)==131 & info(:,4)==65,1),'bo');
plot(pows(info(:,1)==132 & info(:,4)==65,1),'gv');

avPow(2) = mean(pows(info(:,1)==129 & info(:,4)==65,1));
avPow(4) = mean(pows(info(:,1)==130 & info(:,4)==65,1));
avPow(3) = mean(pows(info(:,1)==131 & info(:,4)==65,1));
avPow(1) = mean(pows(info(:,1)==132 & info(:,4)==65,1));

sDev(2) = std(pows(info(:,1)==129 & info(:,4)==65,1));
sDev(4) = std(pows(info(:,1)==130 & info(:,4)==65,1));
sDev(3) = std(pows(info(:,1)==131 & info(:,4)==65,1));
sDev(1)= std(pows(info(:,1)==132 & info(:,4)==65,1));


line([0 45],[avPow(1) avPow(1)],'color','g','LineStyle','-.','LineWidth',1.5);
line([0 45],[avPow(2) avPow(2)],'color','r','LineStyle','-.','LineWidth',1.5);
line([0 45],[avPow(3) avPow(3)],'color','b','LineStyle','-.','LineWidth',1.5);
line([0 45],[avPow(4) avPow(4)],'color','m','LineStyle','-.','LineWidth',1.5);
ylim([-14.5 -12]);

figure;
errorbar(1,avPow(1),sDev(1),'gv');hold on;
errorbar(2,avPow(2),sDev(2),'rx');
errorbar(3,avPow(3),sDev(3),'bo');
errorbar(4,avPow(4),sDev(4),'m^');

xlim([0 5]);
ylim([-14.5 -12]);

ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'Right','Up','Left','Down'};
