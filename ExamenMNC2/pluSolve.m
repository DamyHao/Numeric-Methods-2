function x = pluSolve(L, U, P, b)
    % pluSolve Fa les permutacions i resol el sistema LUx = b
    %   INPUT:
    %       Tot el que treiem de PLU i el vector de termes independents.
    %   OUTPUT:
    %       x = solucio
    n = length(b);

    for i = 1:1:n
        b([i P(i)]) = b([P(i) i]);
    end

    % LUx = b --> Ly = b (trobar y)
    y = FS(L, b);
    % U*x = y
    x = BS(U, y);
end
