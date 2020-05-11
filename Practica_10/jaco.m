function DF = jaco(F,x)
% This code give you the Jacobian matrix of the function F evaluated in x.
% The Jacobian matrix is (n,m) meanwhile the sixe of F is n and the size of
% x is m.

% F son les n funcions escalars


% Input: F(x):R^m ---> R^n
       % x: (m x 1)-vector ; F: (n x 1)-vector
% Output: DF(x) (n x m) Jacobian matrix at x

f1=feval(F,x);    m=length(x);    n=length(f1);

h=sqrt(eps);  H=eye(m)*h;   

DF = zeros(n,m);


for j=1:m
    
    f2=feval(F,x+H(:,j));  
    
    DF(:,j)=(f2-f1)/h;
end



end