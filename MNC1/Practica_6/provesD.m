clear all;
close all
format long g;
z = [0+10^-12:1/2000:1-10^-12];
n=400;
j=[0:1:n];
zj = cos((j*pi)/n);
xj = (1/2)*(1+zj);
n = length(xj);
m = length(z);
den = zeros(length(z),1);
num = zeros(length(z),1);
funxj = fun(xj);
for j = [1:1:n]
        if j == 1 || j == n
        numj = ((-1)^j)*funxj(j)*(z-xj(j)).^(-1).*(1/2);
        num = num + numj';
        denj = ((-1)^j)*(z-xj(j)).^(-1).*(1/2);
        den = den + denj';
    else
        num = num + (((-1)^j)*funxj(j)*(z-xj(j)).^(-1))';
        den = den + (((-1)^j)*(z-xj(j)).^(-1))';
    end      
end

b = num./den;

%b = baricentrica2(z, xj, fun(xj));
figure;
plot(z, b);
figure;
plot(xj, fun(xj));