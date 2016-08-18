function [ band_powers,fBounds ] = create_freqBands( spectrogram, f, bounds )
%create_freqBands - 11/15
% Brendan Thorn
%
% Takes a freq vs. time spectrogram matrix and returns a matrix with
% average power within each band. f, bounds input vectors are optional -
% should be inputted as scalar values if not included. Default frequency
% bands will be returned if vectors are not inputted.
%
% Inputs:
%
%   spectrogram - matrix whose elements are power for given time and frequency
%   f - vector that gives frequency at which power is estimated
%   bounds - vector that defines frequency band upper-bounds for
%   matrix that will be returned.
%
% Outputs:
%
%   band_powers - matrix whose elements are averaged power within given
%   frequency band at each time index
%   fBounds - vector that gives upper frequency-bounds for each frequency
%   band of band_powers matrix
%

% IF no bounds vector given, build band_powers in a default setting
if length(bounds)==1
    % build band_powers matrix with delta(1-4 Hz), theta(4-8 Hz), alpha(8-13 Hz),
    % beta(13-30 Hz), gamma(30-70 Hz), high gamma(70-200 Hz) bands
    band_powers = zeros(2,size(spectrogram,2),size(spectrogram,3));
    
    for i = 1:size(spectrogram,3)
    
        band_powers(1,:,i) = mean(spectrogram(1:7,:,i));
        band_powers(2,:,i) = mean(spectrogram(8:25,:,i));
        %band_powers(3,:,i) = mean(spectrogram(5:7,:,i));
        %band_powers(4,:,i) = mean(spectrogram(8:15,:,i));
        %band_powers(5,:,i) = mean([spectrogram(16:30,:,i); spectrogram(33:36,:,i)]);
        %band_powers(6,:,i) = mean([spectrogram(37:61,:,i); spectrogram(63:92,:,i); spectrogram(94:102,:,i)]);
        
        if i == 1
            fBounds(1) = f(8);
            fBounds(2) = f(26);
            %fBounds(3) = f(8);
            %fBounds(4) = f(16);
            %fBounds(5) = f(37);
            %fBounds(6) = f(103);
        end
    
    end
    
% ELSE bounds are given for frequency bands, so build band_powers
% accordingly
else
    if length(size(spectrogram)) == 3
        % build band_powers matrix with bands specified by bounds vector
        band_powers = zeros(length(bounds),size(spectrogram,2),size(spectrogram,3));
        fBounds = zeros(length(bounds),1);
        % Determine freq bounds based off approximate input bounds
        freqInd = 1;
        for i = 1:length(bounds)
            [~, idx] = min(abs(f-bounds(i)));
            fBounds(i) = f(idx); 
            for j = 1:size(spectrogram,3)
                band_powers(i,:,j) = mean(spectrogram(freqInd:idx-1,:,j));
            end
            freqInd=idx;
        end
    elseif length(size(spectrogram))== 4
        % build band_powers matrix with bands specified by bounds vector
        band_powers = zeros(size(spectrogram,1),length(bounds),size(spectrogram,3),size(spectrogram,4));
        fBounds = zeros(length(bounds),1);
        % Determine freq bounds based off approximate input bounds
        freqInd = 1;
        for i = 1:length(bounds)
            [~, idx] = min(abs(f-bounds(i)));
            fBounds(i) = f(idx); 
            for j = 1:size(spectrogram,4)
                band_powers(:,i,:,j) = mean(spectrogram(:,freqInd:idx-1,:,j),2);
            end
            freqInd=idx;
        end
    end
end

end