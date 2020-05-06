function F = newtonFunction(f)

p = @(x)(x.*0 +1); 
q =@(x)(0.*x + 0);
r = @(x)(0.*x + 0);
n=26;
C=[1 0 0; 1 0 0];

a=0; b=1;
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
mCoef = [ C(2,2), 0 ;  0, C(1,2)];
N=-exp(f+1);
F = (Lhat + M3*inv(M2)*mCoef*M1)*f + M3*inv(M2)*[C(2,3);C(1,3)] - N;

end