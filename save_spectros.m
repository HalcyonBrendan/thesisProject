%%Brendan Thorn - 6/15
% save_spectros.m
% 
% Produces and saves spectrograms of LFP data for each channel. Uses
% MTSpectrogram.m function from FMA Toolbox 
% (see http://fmatoolbox.sourceforge.net/API/FMAToolbox/Analyses/MTSpectrogram.html),
% which also uses the chronux toolbox (see chronux.org).

%save spectrograms from every channel
for j = 1:48
    %construct filename of matfile with channel data in it
    if j < 10
        filename = sprintf('../data/M20100301_246/LFPChan0%d.mat',j);
    else
        filename = sprintf('../data/M20100301_246/LFPChan%d.mat',j);
    end
    
    %open/load matfile into workspace (should overwrite previous one)
    load(filename);
    
    %spectroMag = zeros(257,578,590);
    spects = zeros(257,369,590);
    t = zeros(590,369);
    f = zeros(590,257);
    data = zeros(length(LFP.AD(1,:)),2);

    %produce and save spectrogram matrices
    for i = 1:590
        %fft spectrograms
        %spectroMag(:,:,i) = spectrogram(LFP.AD(i,:),256,240,512,1e3,'yaxis');

        %multitaper spectrograms
        data(:,1) = double(LFP.AD(i,:));
        data(:,2) = double(LFP.AD(i,:));
        
        %uses NW=3, K=5 by default - TODO: look at NW=4, K=5 among others
        [spects(:,:,i),t(i,:),f(i,:)] = MTSpectrogram(data,'frequency',1000,'window',.3,'step',.025);
        
    end

    handle = sprintf('../analysis/mt_spectrograms/M20100301_246/mt_spects_ch%d.mat',j);

    save(handle,'spects','t','f');

    
    disp(j);
    
end