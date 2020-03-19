function d = der(f,Xo);
h= 10.^-10;
f1 = f(Xo+h);
f2 = f(Xo);
d = (f1-f2)/h;
end


