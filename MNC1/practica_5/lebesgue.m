function l = lebesgue(n)
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    z = [-1:1/300:1];
    P = MLagrange(z, x);
    for j=0:n
        l = sum(abs(P),2); %suma les files de la matriu P en cada cas.
    end
end
