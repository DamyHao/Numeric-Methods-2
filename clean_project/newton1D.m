function xkfin = newton1D(a, b, tol, itmax, fun)
it=0;
tolk=1;
xk = [(a+b)/2];

if fun(a)*fun(b)<0
    while (it<itmax && tolk>tol)
        xk = [xk xk(end)-(fun(xk(end))/der(fun,xk(end)))];
        tolk = abs(xk(end-1)-xk(end));
        it = it+1;
    end
%else
    %disp('les imatges de a i b tenen signes iguals');
end
xkfin = xk(end);
end

%%
function d = der(fun, Xo)
h = 10^-10;
d = (fun(Xo+h)-fun(Xo))/h;
end
          