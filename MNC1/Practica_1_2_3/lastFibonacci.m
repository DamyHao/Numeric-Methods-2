function out = lastFibonacci(ultim)
if ultim <= 0
    out = 0;
elseif ultim == 1 || ultim == 2
    out = 1;
else
    n1 = 1;
    n2 = 1;
    temp = 2;
    aux = 0;
    while temp < ultim
        aux = n1 + n2;
        n1 = n2;
        n2 = aux;
        temp = temp + 1;
    end
    out = n2;
end
    
    