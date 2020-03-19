function out = fibonacci(n)
if n <= 0
    out = 0;
elseif n == 1 || n == 2
    out = 1;
else
    out = fibonacci(n-1) + fibonacci(n-2);
end
