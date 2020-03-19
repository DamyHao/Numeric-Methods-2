function [xk,err,it,p] = newton(a,b, tol, itmax, fun)
it=0; 
iterr = 1;
tolk=1;

xk = [(a+b)/2];

syms x;
funder = diff(sym(fun), 1);

    if fun(a)*fun(b)<0
        while (it<itmax && tolk>tol)
            xk = [xk xk(end)-(fun(xk(end))/der(fun,xk(end)))];
            tolk = abs(xk(end-1)-xk(end));
            it = it+1;
        end
        err = [];
        while iterr < it + 2
            err = [err abs(xk(iterr) - xk(end))];
            iterr = iterr + 1;
        end
                
    else
        disp('les imatges de a i b tenen signes iguals');
    end
    p = 0;
    
end