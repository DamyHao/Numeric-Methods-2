% Function that solves linear system of the form Ax=b with A factorized as PA=LU
% Input:  L lower triangular nxn matrix
%         U upper triangular nxn matrix
%         P reoredering vector nx1 Rn vector
%         b right hand side of the system nx1 Rn vector
%
% Output: solution x
% Note:   requires fs.m and bs.m
function x = solvePPLU(L,U,P,b)
    n=length(b);
    for k = 1:n-1
        b([k P(k)]) = b([P(k) k]);
    end
    y = fs(L,b);
    x = bs(U,y);
end