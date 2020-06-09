% Function that performs LU factorization with partial pivoting on a matrix A, LU=PA
% Input:    A nxn matrix
%
% Output:   L lower triangular nxn matrix
%           U upper triangular nxn matrix
%           P reoredering vector nx1 Rn vector
function [P,L,U]=factPPLU(A)
    n=length(A);
    U=A; L=eye(n);P=[1:n]';
    for i=1:n-1
        [maxvalue,maxpos]=max(abs(U(i:n,i)));
        maxpos=maxpos+i-1;
        i1=[i,maxpos];
        i2=[maxpos, i];
        U(i1,:)=U(i2,:);
        P(i)=maxpos;
        L(i1,1:i-1) = L(i2,1:i-1);
        for k=i+1:n
            L(k,i)=U(k,i)/U(i,i);
            U(k,i)=0;
            for j=i+1:n                         
                U(k,j)=U(k,j)-L(k,i)*U(i,j);
            end
        end
    end
end