function [H]=hessV(X)
x=X(1,1);
y=X(2,1);

h1=[ 12*x^2 + 4*y - 42,         4*x + 4*y];
h2=[         4*x + 4*y, 12*y^2 + 4*x - 26];
H=[h1;h2];
end