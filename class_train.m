
function [ NBModel ] = class_train(X,y)
%class_train - 11/15
% Brendan Thorn
%
% Train a naive bayesian classifier with input:
%
% X - size(X)=[n,m], n training samples, each of which consists of 
% m=numTimeSteps*numFreqBand data points
% y - siz(y)=[n,1], class labels for each training sample, which ultimately
% specify a unique combination of movement type and direction

NBModel = fitNaiveBayes(X,y,'Prior','uniform');

end