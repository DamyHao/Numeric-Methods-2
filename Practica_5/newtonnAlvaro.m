% Code 20: Newton’s method for n-dimensional systems
% Input: x0 - initial guess (column vector)
% tol - tolerance so that ||x_{k+1} - x_{k} || < tol
% itmax - max number of iterations
% fun - function’s name
% Output: XK - iterated
% resd: resulting residuals of iteration: ||F_k||
% it: number of required iterations to satisfy tolerance
function [XK, resd, it] = newtonnAlvaro(x0, tol, itmax, fun)
    xk = [x0]; resd = [norm(feval(fun, xk))]; XK = [x0]; it = 1;
    tolk = 1.0; n = length(x0);

    while it < itmax & tolk > tol
        Fk = feval(fun, xk);
        DFk = jaco(fun, xk); [P, L, U] = PLU(DFk);
        %dxk = pluSolve(L, U, P, -Fk);
        dxk = DFk\(-Fk)';
        xk = xk + dxk; XK = [XK xk]; resd = [resd norm(Fk)];
        tolk = norm(XK(:, end) - XK(:, end - 1)); it = it + 1;
    end
