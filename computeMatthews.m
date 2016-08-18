function [ mcc ] = computeMatthews( tp,tn,fp,fn )
%computeMatthews.m - 03/16 Brendan Thorn
%   computeMatthews.m computes the Matthews correlation coefficient given
%   number of true positives tp, true negatives tn, false positives fp,
%   false negatives fn, for a binary classifier.

mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));

if isnan(mcc)
    mcc = 0;
elseif isinf(mcc)
    mcc = 0;
end

end