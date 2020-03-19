function [xk,it,p] = newtonD(a,b, tol, itmax, fun)
it=0; 
tolk=1;
xk = [(a+b)/2];

    if fun(a)*fun(b)<0
        while (it<itmax && tolk>tol)
            xk = [xk xk(end)-(fun(xk(end))/derD(fun,xk(end)))];
            tolk = abs(xk(end-1)-xk(end));
            it = it+1;
        end
    else
        disp('les imatges de a i b tenen signes iguals');
    end
    p = 0;
end