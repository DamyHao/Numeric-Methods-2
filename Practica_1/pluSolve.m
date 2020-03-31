

function x = pluSolve(L, U, P, b)
    % pluSolve Fa les permutacions i resol els sistemes
    n = length(b);

    for i = 1:1:n
        b([i P(i)]) = b([P(i) i]);
    end

    y = FS(L, b);
    x = BS(U, y);

end
