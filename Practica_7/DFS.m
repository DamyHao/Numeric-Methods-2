function FNx = DFS (fk,x)
% Evaluada en un malla de punts x la serie de fourier a partir dels
% coeficients de fourier fk. x HA DE SER una columna

N = length(fk); M = length(x) ; k = [0:N-1]; 
PHI = exp(i*x*(k-N/2)); FNx = PHI*fk;

end