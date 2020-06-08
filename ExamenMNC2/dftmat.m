function fk = dftmat(fj)
    % Li donem els punts avaluats en physical space i ens retorna els
    % coeficients de la serie de Fourier (transformades de fourier (fourier space))
    % tot en columna
%     N = length(fj);
%     W = exp(-1i * 2 * pi / N);
%     k = 0:N - 1; j = 0:N - 1;
% 
%     F = (1 / N) * W.^((k' - N / 2) * j);
%     fk = F * fj;

N = length(fj); WN = exp(-1i*2*pi/N); jj = [0:N-1]; kk = jj';
F = (1/N)*WN.^((kk-N/2)*jj); fk = F*fj;
end
