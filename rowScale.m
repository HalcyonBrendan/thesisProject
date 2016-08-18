function [ Xout ] = rowScale( Xin,dim2 )
%rowScale - May 3,2016 Brendan Thorn
%   Accepts matrix Xin and multiplicatively scales/normalizes the rows such
%   that the sum of the absolute value of each row equals 1000. Only does
%   so if the dimensionality of Xin is correct though, which cannot
%   necessarily be determined just from Xin (for thesis applications).
    
    Xout = Xin;
    if dim2 > 1
        Xsums = sum(abs(Xin),2);
        factors = 1000./Xsums;
        Xout = Xin.*repmat(factors,1,size(Xin,2));
    end

end

