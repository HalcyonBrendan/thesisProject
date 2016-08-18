function [ yOut,score ] = svm_test(SVMModel,Xtest)
%svm_test - 02/16
% Brendan Thorn
%
% Classifies data given in Xtest according to SVM model trained in
% svm_train.m and given by SVMModel data structure

[yOut,score] = predict(SVMModel,Xtest);

end