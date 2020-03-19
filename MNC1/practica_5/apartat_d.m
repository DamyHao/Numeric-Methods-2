% intent apartat d)
clear all;
z = [-1:1/300:1];
    n=4; % nomes provisional%
i_exp = interpolacio_exp(n);
figure(10);
plot(z,i_exp,'<r');
hold on;
