% Code 22: DFT (Matrix-vector product)
% Input: fj (N-column vector: sampling values of f)
% Output: fk (N-column vector: Fourier coefficients)
function fk = dftmat(fj)
N = length(fj); WN = exp(-i*2*pi/N); jj = [0:N-1]; kk = jj';
F = (1/N)*WN.^((kk-N/2)*jj); fk = F*fj;
end