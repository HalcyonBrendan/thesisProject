%%Brendan Thorn - 12/15
% mse_classification.m
% 
% Performs several classification tasks, such as determining directional
% tuning during eye and reach movements, determining behavioural state (ie.
% reaching, eye movement, planning). For BT M.Eng Thesis.
% 


% Use minimum mean square error during mem period to classify eye
% movements into a category defined by type (saccade or pursuit)
% and direction (right, up, left, down).

fprintf('Beginning classification2.m!\n');

% Load an LFP data structure into workspace for trial info
filename = sprintf('../data/monkey_data/LFPChan01.mat');
load(filename);
% For later convenience and code clarity, give data summary matrix a shorter name
info = Events.dataSummary;

% Get total number of trials from given data
numTrials = size(info,1);

% Choose channels that you wish to use for classification
channels = 1:48;
numChannels = length(channels);

% How many time steps we want to look at for each band
numTimeSteps = 3;

% Before starting, a bit of code to put eye move times into one vector
alignTimes = zeros(numTrials,1);
for i = 1:numTrials
    if info(i,4)==65
        alignTimes(i)=TT.T65sac(i);
    elseif info(i,4)==66
        alignTimes(i)=TT.T66on(i);
    else
        alignTimes(i)=TT.T8(i);
    end
end

