function D = diferenciacion(x)
llarg = length(x);
D = zeros(llarg, llarg);
for i = 1:1:llarg
    for j = 1:1:llarg
        if i ~= j
            D(i,j) = (lambda(x, j)/lambda(x, i)) * 1/(x(i) - x(j));
        else
            sum = 0;
            for k = 1:1:llarg
                if k ~= j
                    sum = sum + (1/(x(j)-x(k)));
                end
            end
            D(i,j) = sum;
        end
    end
end
end