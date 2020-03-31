function [y,iconv] = continuationStep(fun,y0,y1,s,tol,itmax)


    resd = [norm(feval(fun, xk))]; 
    it = 1; 
    tolk = 1;
    v=y1-y0;
    xk = y1 + v*s % Si s = 1 conseguim que la separació entre solucions sigui el maxim de "constant"

    % A part de les ecuacions que teniem en el nnewton normal, li
    % imposareem que el preoducte escalar entre v i (xk(punt
    % buscat)- xk(predictor)) sigui 0
    while it < itmax && tolk > tol
        J = jaco(fun, xk); % Jacobia en la posicio anterior
        J = [J ; v];
        fk =[ -fun(xk)'; v*' ]; 
        %[P, L, U] = PLUAlvaro(J);
       
        %Dx = pluSolve(L, U, P, (-fk)'); %Solucio de la ecuacio J*Dx = -fk

        Dx = J\fk;
        xk = xk + Dx;
        XK = [XK, xk];
        resd = [resd, norm(fk)];
        tolk = norm(XK(:, end) - XK(:, end - 1));
        it = it + 1;
        

    end




end