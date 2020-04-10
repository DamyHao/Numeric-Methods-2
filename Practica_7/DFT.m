function fk=DFT(fj)
% Li donem els punts avaluats en physical space i ens retorna les
% transformades de fourier (fourier space)
%fj que sigui columna
N=length(fj);
W=exp(-i*2*pi/N);
k=0:N-1; j=0:N-1;

F=(1/N)*W.^((k'-N/2)*j);
fk = F*fj;

end