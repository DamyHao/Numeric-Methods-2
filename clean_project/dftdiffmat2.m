% Code 26: Fourier 2nd Ord Differentiation Matrix % Input: N (even)
% Output: NxN Diff. Matrix D2
function D2 = dftdiffmat2(N)
% Segona derivada d'una funcio 2pi periodica.
% Cal multiplicar aquesta matriu per una malla de la funcio que volem
% diferenciar avalauda en els nodes xj = 2pi/N.
WN = exp(-1i*2*pi/N) ; D2 = zeros(N) ; k = [0:N-1];
for j = 0:N-1
    for l = 0:N-1
        D2(j+1,l+1) = -sum(((k-N/2).^2).*WN.^(-(k-N/2)*(j-l)))/N;
    end
end