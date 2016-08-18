function [ normBands ] = normByBand( spectrogram,band )
%normByBand.m - 010216
%   NormByBand takes a freqxtime spectrogram and normalizes each column by
%   dividing by the value of the element in that column and in the row
%   specified by band.

normalizer = spectrogram(band,:,:);
normBands = spectrogram./repmat(normalizer,size(spectrogram,1),1,1);

end

