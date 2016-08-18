mip6 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip6_dirSac_p01_250216.mat');
mip6mpc = mip6.meanPercentCorr;
mip6numTrainTrials = mip6.trainNums;

x=mip6numTrainTrials';
y=mip6mpc';

myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);

figure;
p=plot(myfit,'g',x,y,'g*');
title('Channel Comparison');
hold on;

mip12 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_p01_250216.mat');
mip12mpc = mip12.meanPercentCorr;
mip12numTrainTrials = mip12.trainNums;

x=mip12numTrainTrials';
y=mip12mpc';

myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);

p=plot(myfit,'k',x,y,'kx');
title('Channel Comparison');
hold on;


mip12_0 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_0.mat');
mip12_1 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216.mat');
mip12_2 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_2.mat');
mip12_3 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_3.mat');
mip12_4 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_4.mat');
mip12_5 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_5.mat');
mip12_6 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_6.mat');
mip12_7 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_7.mat');
mip12_8 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_250216_8.mat');
mip12mpc = [];
mip12mpc = [mip12mpc mip12_0.meanPercentCorr];
mip12mpc = [mip12mpc mip12_1.meanPercentCorr];
mip12mpc = [mip12mpc mip12_2.meanPercentCorr];
mip12mpc = [mip12mpc mip12_3.meanPercentCorr];
mip12mpc = [mip12mpc mip12_4.meanPercentCorr];
mip12mpc = [mip12mpc mip12_5.meanPercentCorr];
mip12mpc = [mip12mpc mip12_6.meanPercentCorr];
mip12mpc = [mip12mpc mip12_7.meanPercentCorr];
mip12mpc = [mip12mpc mip12_8.meanPercentCorr];
mip12numTrainTrials = [];
mip12numTrainTrials = [mip12numTrainTrials mip12_0.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_1.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_2.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_3.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_4.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_5.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_6.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_7.trainNums];
mip12numTrainTrials = [mip12numTrainTrials mip12_8.trainNums];

x=mip12numTrainTrials';
y=mip12mpc';

myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);

p=plot(myfit,'b',x,y,'ob');
hold on;


mip18_0 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_0.mat');
mip18_1 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216.mat');
mip18_2 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_2.mat');
mip18_3 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_3.mat');
mip18_4 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_4.mat');
mip18_5 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_5.mat');
mip18_6 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_6.mat');
mip18_7 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_7.mat');
mip18_8 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_250216_8.mat');
mip18mpc = [];
mip18mpc = [mip18mpc mip18_0.meanPercentCorr];
mip18mpc = [mip18mpc mip18_1.meanPercentCorr];
mip18mpc = [mip18mpc mip18_2.meanPercentCorr];
mip18mpc = [mip18mpc mip18_3.meanPercentCorr];
mip18mpc = [mip18mpc mip18_4.meanPercentCorr];
mip18mpc = [mip18mpc mip18_5.meanPercentCorr];
mip18mpc = [mip18mpc mip18_6.meanPercentCorr];
mip18mpc = [mip18mpc mip18_7.meanPercentCorr];
mip18mpc = [mip18mpc mip18_8.meanPercentCorr];
mip18numTrainTrials = [];
mip18numTrainTrials = [mip18numTrainTrials mip18_0.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_1.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_2.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_3.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_4.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_5.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_6.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_7.trainNums];
mip18numTrainTrials = [mip18numTrainTrials mip18_8.trainNums];

x=mip18numTrainTrials';
y=mip18mpc';

myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);

plot(myfit,'r',x,y,'dr');
hold on;

mip36 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip36_dirSac_p01_250216.mat');
mip36mpc = mip36.meanPercentCorr;
mip36numTrainTrials = mip36.trainNums;

x=mip36numTrainTrials';
y=mip36mpc';

myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);

plot(myfit,'c',x,y,'c+'); xlim([0 175]);ylim([.2 .5]);
xlabel('Number of Training Trials');ylabel('Accuracy [%]');
hold off;
%%%%%%%%%%%
%%
%%%%%%%%%%%
mip2 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip2_dirSac_010316.mat');
mip2mpc = mip2.meanPercentCorr;
mip2numTrainTrials = mip2.trainNums;
x=mip2numTrainTrials';
y=mip2mpc';
myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);
figure;
p=plot(myfit,'g',x,y,'g*');
hold on;

mip4 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip4_dirSac_010316.mat');
mip4mpc = mip4.meanPercentCorr;
mip4numTrainTrials = mip4.trainNums;
x=mip4numTrainTrials';
y=mip4mpc';
myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);
p=plot(myfit,'b',x,y,'b+');
hold on;

mip8 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip8_dirSac_010316.mat');
mip8mpc = mip8.meanPercentCorr;
mip8numTrainTrials = mip8.trainNums;
x=mip8numTrainTrials';
y=mip8mpc';
myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);
p=plot(myfit,'r',x,y,'rx');
hold on;

mip12 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip12_dirSac_010316.mat');
mip12mpc = mip12.meanPercentCorr;
mip12numTrainTrials = mip12.trainNums;
x=mip12numTrainTrials';
y=mip12mpc';
myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);
p=plot(myfit,'k',x,y,'kv');
hold on;

mip18 = load('/Users/brendan/Documents/Master/Project/LFPs/analysis/workspaces/M20100302_645/w_mip18_dirSac_010316.mat');
mip18mpc = mip18.meanPercentCorr;
mip18numTrainTrials = mip18.trainNums;
x=mip18numTrainTrials';
y=mip18mpc';
myfittype = fittype('a + b*log(x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});
myfit = fit(x,y,myfittype);
p=plot(myfit,'c',x,y,'cs');
title('Channel Comparison');
hold off;
xlabel('Number of Training Trials'); ylabel('Decode Accuracy [%]');
xlim([0 175]);ylim([.22 .55]);
