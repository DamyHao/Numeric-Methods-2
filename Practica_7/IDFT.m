function fj= IDFT (fk)
% Et dona la funcio en el physicial space evaluada en els punts xj a partir
% dels coeficients de fourier (en format columna)

N = length(fk); WN = exp(-i*2*pi/N); jj = [0:N-1]'; kk = jj'; 
P = WN.^(-jj*(kk-N/2)); fj = P*fk

end