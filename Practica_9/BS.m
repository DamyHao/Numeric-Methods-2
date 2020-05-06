function x = BS(A, b)
    n = length(A);
    x = zeros(n, 1);

    for i = n:-1:1
        sum = 0;

        for j = i + 1:n
            sum = sum + A(i, j) * x(j);
        end

    end

end
