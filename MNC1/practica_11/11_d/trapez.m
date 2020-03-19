function I_trap = trapez(a, b, fun)
%TRAPEZ Summary of this function goes here
%   Detailed explanation goes here

h=b-a;
I_trap = (h/2)*(fun(a)+fun(b));
end

