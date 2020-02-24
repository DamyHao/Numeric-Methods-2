function x = BS(UT, b) 

[dim, n] = size(UT);
if dim~= n; error('NECESSITA UNA MATRIU QUADRADA. Revisar els arguments'); end

x=[0*(1:n-1),b(n)/UT(n,n)];

for i=n:-1:1
    sum = 0;
    
    for j=n:-1:i+1
        sum = sum + UT(i,j)*x(j);  
    end
    x(i)=(1/UT(i,i))*(b(i)-sum);
end
end

