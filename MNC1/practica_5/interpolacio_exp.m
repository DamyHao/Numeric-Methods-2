function [exp_x,i_exp] = interpolacio_exp(n)
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    z = [-1:1/300:1];
    P = MLagrange(z, x);
    numero_z = length(z); %dubto si ha de ser length z o x
    for j=1:numero_z
    exp_x = exp(x);
    P_sumatori = P(j,:);
    prod = exp_x.*P_sumatori;
    i_exp = sum(prod);
    end
end

    