function [ SVMModel ] = svm_train(X,y)
%svm_train - 02/16
% Brendan Thorn
%
% Train a SVM classifier:

t = templateSVM('Standardize',1,'KernelFunction','polynomial','PolynomialOrder',4,'KernelScale','auto');
SVMModel = fitcecoc(X,y,'Learners',t,'Prior',[.25 .25 .25 .25],'FitPosterior',1);

end