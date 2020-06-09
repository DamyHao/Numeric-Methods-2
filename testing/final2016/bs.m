% Function that performs backwards substitution to solve system Ax=b
% Input:    A upper triangular nxn matrix
%           b R^n column vector
%
% Output:   x R^n column vector  
function x=bs(A,b)
    n=length(A);
    x=zeros(n,1);
    for i=n:-1:1
        sum=0;
        for j=i+1:n
            sum=sum+A(i,j)*x(j);
        end
        x(i)=(1/A(i,i))*(b(i)-sum);
    end
end