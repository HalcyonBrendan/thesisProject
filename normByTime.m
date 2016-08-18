function [ normBands ] = normByTime( freqBands )
%normByTime.m - 02/16 Brendan Thorn
%   normByTime.m takes frequency band-averaged log-spectrograms and returns
%   normBands, which is freqBands scaled such that the sum of the average
%   log-powers in the bands at each given time is equal to -100.

powerSum = sum(freqBands,1);
normBands = 100*freqBands./repmat(powerSum,size(freqBands,1),1);

end

