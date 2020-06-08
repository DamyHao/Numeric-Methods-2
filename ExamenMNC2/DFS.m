function FNx = DFS (fk,x)
% Evaluada en un malla de punts x la serie de fourier a partir dels
% coeficients de fourier fk. (simplement multiplica els coefiecients per lespai de funcions)
% x HA DE SER una columna

N = length(fk); k = [0:N-1]; 
PHI = exp(1i*x*(k-N/2)); FNx = PHI*fk; 
% PHI sera matriu:
% [ phi_-N/2(xo)     ...     phi_N/2-1(xo) ]
% [   ...            ...            ...    ]
% [ phi_-N/2(x_n-1)  ...  phi_N/2-1(x_n-1) ]
% On cada phi es un un vector de la base unitaria de l'espai de funcions
% phi_a(x) = e^(i*a*x)

end