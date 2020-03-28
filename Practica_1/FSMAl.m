function x = FS(LT, b)
format long
[dim, n] = size(LT);
if dim~= n; error('NECESSITA UNA MATRIU QUADRADA. Revisar els arguments'); end

x=[b(1)/LT(1,1),0*(1:n-1)];

for i=1:1:n
    sum = 0;
    
    for j=1:1:i-1
        sum = sum + LT(i,j)*x(j);  
    end
    x(i)=(1/LT(i,i))*(b(i)-sum);
end

end