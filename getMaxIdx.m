% getMaxIdx.m
% Brendan Thorn - 210316
%
% Returns the maximum value and the index vector of that value for a
% multidimensional matrix. Right now I am just going to do a simple, ugly
% version for 4D matrices. I'll hopefully do a better, more general job
% later.

function [ val, idx ] = getMaxIdx( mat )

% Get size and number of dimensions of mat
sizes = size(mat);
% Keep track of the results
val = 0; idx = zeros(3,1);

for i = 1:sizes(1)
    for j = 1:sizes(2)
        for k = 1:sizes(3)
            if mat(i,j,k) > val
                val = mat(i,j,k);
                idx(1)=i;idx(2)=j;idx(3)=k;
            end
        end
    end
end

end