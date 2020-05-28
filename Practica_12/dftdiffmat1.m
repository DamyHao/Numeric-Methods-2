% Code 25: Fourier 1st Ord Differentiation Matrix % Input: N (even)
% Output: NxN Diff. Matrix D1 %ALVARO
function D1 = dftdiffmat1(N)
% Deriva una funcio 2pi periodica
% Cal multiplicar aquesta matriu per una malla de la funcio que volem
% diferenciar avalauda en els nodes xj = 2pi/N
WN = exp(-i*2*pi/N) ; D1 = zeros(N) ; k = [0:N-1]; 
for j = 0:N-1
    for l = 0:N-1
        D1(j+1,l+1) = (i/N)*sum((k-N/2).*WN.^(-(k-N/2)*(j-l)));
    end
end