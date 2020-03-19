function xkfin = newton2(a, tol, itmax, fun)
%NEWTON2 Summary of this function goes here
%   Newton que nomes es llança en un punt.

it=0;
tolk=1;
xk = [a];

    while (it<itmax && tolk>tol)
        xk = [xk xk(end)-(fun(xk(end))/der(fun,xk(end)))];
        tolk = abs(xk(end-1)-xk(end));
        it = it+1;
    end

xkfin = xk(end);
end

