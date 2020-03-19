function d = der(Xo,f)
h= 10^-10;
d = (f(Xo+h)-f(Xo))/h;
end


