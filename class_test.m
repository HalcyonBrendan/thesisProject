function [ yOut ] = class_test(NBModel,Xtest)
%class_test - 11/15
% Brendan Thorn
%
% Classifies data given in Xtest according to naive bayes model trained in
% class_train.m and given by NBModel data structure

yOut = predict(NBModel,Xtest);

end