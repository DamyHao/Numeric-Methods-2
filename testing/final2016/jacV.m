function [Y]=jacV(X)
x=X(1,1);
y=X(2,1);
Y=[ 2*x + 4*x*(x^2 + y - 11) + 2*y^2 - 14; 2*y + 4*y*(y^2 + x - 7) + 2*x^2 - 22];
end