function [ normSpect ] = normalize_spect(spectrogram,idxs)
%normalize_spect - 11/15
% Brendan Thorn
%
% Normalize a spectrogram by its baseline activity with input:
%
% spectrogram - (3d) matrix to be normalized
% idxs - numTrialsxm matrix with indices for each trial that comprise 
% baseline activity in spectrogram and will be averaged to provide 
% baseline power for each trial.

normSpect = zeros(size(spectrogram,1),size(spectrogram,2),size(spectrogram,3));

for i = 1:size(spectrogram,3)
    
    baselineMean = mean(spectrogram(:,idxs(i,:),i),2);
    
    normSpect(:,:,i) = s(spectrogram(:,:,i))./repmat(baselineMean,1,size(spectrogram,2));

end