%function [  ] = spectrogram_bt1( x, timeStep, windowLength )
% spectrogram_bt1.m
%
% Brendan Thorn
% May 21, 2015
%
% Function to produce spectrograms (heat maps of freq vs. time) with
% multitaper method for power spectrum estimation

%delete this block if using as a function
%%%%%
x = LFP.AD(1,:);
timeStep = 50;
windowLength = 300;
%%%%%

length = windowLength;
step = timeStep;

nfft = 2^nextpow2(length);
nw=4; %default
fs = 1000; %ie. millisecond resolution

here = 1;

numSteps = ceil(size(x,2)/step);

spectro = zeros(257,numSteps);

for stepNum = 1:numSteps-60
    
    spectro(:,stepNum) = pmtm(x(here:here+length-1),nw,nfft,fs);
    
    here = here + step;
    
end

spectrogram = log(spectro);
spec=imagesc(spectrogram);
set(gca,'YDir','normal')
colormap jet;
%end

