function [D,x] = chebdiff(n,a,b)
z =cos([0:n]'*pi/n); d = [.5 ;ones(n-1,1);.5]; % Va de -1 a 1
D = zeros(n+1,n+1);
for ii = 0:n
    for jj = 0:n
        ir = ii + 1 ; jc = jj + 1;
        if ii == jj
            kk = [0:ii-1 ii+1:n]'; num = (-1).^kk.*d(kk+1) ;
            D(ir,jc) =((-1)^(ir)/d(ir))*sum(num./(z(ir)-z(kk+1)));
        else
            D(ir,jc) = d(jc)*(-1)^(ii+jj)/((z(ir)-z(jc))*d(ir));
        end
    end
end
D=2/(b-a) .* D;
x = a+(b-a)/2 * (z+1);% x son nodes Chebyshev amb 'allargats' a el domini del problema