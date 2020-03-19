function d = der(fun, Xo)
h = 10^-10;
d = (fun(Xo+h)-fun(Xo))/h;
end
