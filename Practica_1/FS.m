function x = FS(A, b)
    n = length(A);
    x = zeros(n, 1);

    for i = 1:n
        sum = 0;

        for j = 1:i-1
            sum = sum + A(i, j) * x(j);
        end
        x(i) = (1/A(i,i))*(b(i)-sum);
    end

end
