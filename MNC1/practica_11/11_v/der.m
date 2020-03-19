function d = der(Xo,fun)
h = 10^-10;
d = (fun(Xo+h)-fun(Xo))/h;
end
