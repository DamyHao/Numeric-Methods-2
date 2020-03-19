function x = newton(a, fun)
itmax=20;
tol=10^-6;
tolk=1;
xk=[a a/2];
it=0;
while tolk>tol && it<itmax
    xk= [xk xk(end)-(fun(xk(end))/der(xk(end),fun))];
    it=it+1;
    tolk=abs(xk(end)-xk(end-1));
end
x = xk(end);
end


          