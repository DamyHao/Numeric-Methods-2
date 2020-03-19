function [xk,it,p] = newton(a,b, tol, itmax, fun)
it=0; 
tolk=1;
xk = [(a+b)/2];
   while (it<itmax && tolk>tol)
            xk = [xk xk(end)-(fun(xk(end))/derD(fun,xk(end)))];
            tolk = abs(xk(end-1)-xk(end));
            it = it+1;
   end
   
end