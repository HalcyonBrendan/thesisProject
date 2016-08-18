function [ gram, time, freq ] = getCohs( gramNum, path, chan1, chan2, folder )
%getCohs.m - 11/15 Brendan Thorn
%   getCohs.m uses code from coherogram.m to return one of a spectrogram,
%   coherogram, or cross-power spectrogram. 
%
%   Inputs:
%   1) gramNum - number to specify which gram to return: 1=S1,2=S2,3=S12,4=C
%   2) path - if path is included, return specified gram from same file.
%   Otherwise, path=0 means run code to produce grams from scratch.
%
%   Note: Values for spectrogram generation (eg. number of tapers) are
%   currently set below, but could be taken as inputs in later versions.
%

% If no path is included, assume that the desired matrix must be generated
if length(path) < 2
    % Create arrays of channels of interest (old variable names can be
    % ignored)
    mipChans = chan1;     %1:32;
    pmdChans = chan2;     %33:48

    % Set parameters for cohgramc.m computations later
    % Window parameters [length[s] steplength[s]] 
    % For spect_dec0215_ch%d_tpr%d.mat type spectrograms:
    movingwin = [.3 .025];
    % For spect_mar1716_ch%d_tpr%d.mat type spectrograms:
    %movingwin = [.24 .024];
    % Struct with remaining parameters: tapers, pad, Fs, fpass, err, trialave
    tapers = [4 7]; pad = 1; Fs = 1000; fpass = [0 Fs/2]; err = [0 1]; trialave = 0;
    params = struct('tapers',tapers,'pad',pad,'Fs',Fs,'fpass',fpass,'err',err,'trialave',trialave);

    % For channels of interest in MIP (Chans 1-32 in current data (10/15))
    for j = 1:length(mipChans)
        % Construct filename of matfile with MIP channel data in it
        if mipChans(j) < 10
            filename = sprintf('../data/%s/LFPChan0%d.mat',folder,mipChans(j));
        else
            filename = sprintf('../data/%s/LFPChan%d.mat',folder,mipChans(j));
        end
        % Load matfile into workspace
        mipChan = load(filename);

        % For channels of interest in PMd (Chans 33-48 in current data (10/15))
        for i = 1:length(pmdChans)

            % Construct filename of matfile with PMd channel data in it
            if pmdChans(j) < 10
                filename = sprintf('../data/%s/LFPChan0%d.mat',folder,pmdChans(j));
            else
                filename = sprintf('../data/%s/LFPChan%d.mat',folder,pmdChans(j));
            end
            % Load matfile into workspace
            pmdChan = load(filename);

            % Create matrices that stores LFP data to be processed
            % TO SET: This is only valid for 170316 spects
            data1 = mipChan.LFP.alignedLFP(:,11:end)';
            data2 = pmdChan.LFP.alignedLFP(:,11:end)';

            [C,~,S12,S1,S2,t,f] = cohgramc(data1,data2,movingwin,params);

            %handle = sprintf('../analysis/coherograms/%s/cohgram_nov1015_chs%d_%d_tpr%d.mat',folder,mipChans(j),pmdChans(i),47);
            %save(handle,'C','phi','S12','S1','S2','t','f','confC','phistd','Cerr','-v7.3');

            % Find the correct values to return
            if gramNum == 1
                gram = S1;
            elseif gramNum == 2
                gram = S2;
            elseif gramNum == 3
                gram = S12;
            elseif gramNum == 4
                gram = C;
            else
                gram = 'Error: Could not find desired spectrogram at given path';
            end

            time=t;
            freq=f;
            
            handle = sprintf('../analysis/coherograms/%s/spect_dec0215_ch%d_tpr%d.mat',folder,mipChans(j),47);
            save(handle,'S1','t','f','-v7.3');
        end
    end

% Else if a path is included, load the file and return the desired matrix
else
    % Load file given by path
    load(path);
    % Find the correct matrix
    if gramNum == 1
        gram = S1;
    elseif gramNum == 2
        gram = S2;
    elseif gramNum == 3
        gram = S12;
    elseif gramNum == 4
        gram = C;
    else
        gram = 'Error: Could not find desired spectrogram at given path';
    end

    time=t;
    freq=f;
end

end