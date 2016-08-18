function [ Xstdzed,meanX,stdX ] = stdze( X,dim,meanX,stdX )
%stdze.m - Brendan Thorn, 04/16
% Takes 2D matrix X as input. Standardizes along dimension dim by
% subtracting mean of each row/col from each value, then divides by
% standard deviation of that row/col

if nargin < 3

    meanX = mean(X,dim);
    stdX = std(X,0,dim);

    if dim == 1
        Xstdzed = (X-repmat(meanX,size(X,1),1))./repmat(stdX,size(X,1),1);
    elseif dim == 2
        Xstdzed = (X-repmat(meanX,1,size(X,2)))./repmat(stdX,1,size(X,2));
    end

elseif nargin == 4
    if dim == 1
        Xstdzed = (X-repmat(meanX,size(X,1),1))./repmat(stdX,size(X,1),1);
    elseif dim == 2
        Xstdzed = (X-repmat(meanX,1,size(X,2)))./repmat(stdX,1,size(X,2));
    end
else
    fprintf('Wrong number of inputs to stdze()!');
end

end

