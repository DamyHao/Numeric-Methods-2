function l = lambda(x, j)
% La j es el index fet servir en teoria!!
llarg = length(x);
mult = 1;
for k = 1:1:llarg
    if j ~= k
        mult = mult * (x(j) - x(k));
    end
end

l = 1/mult;
end