%% PRACTICA_5:NEWTON'S METHOD IN R^n
% Damia Casas & Núria Mercadé

%% A) 
% For the temperatures T = {0.99, 0.98, 0.97, . . . , 0.85}, use Newton?s
% method to compute the coordinates vl(T) and vg(T), along with the
% corresponding pressure. When changing T, use the previous result as
% initial guess. For T = 0.99, use v(0) = 0.8 and v(0) = 1.2.

t=[0.99:-0.01:0.85]; 

v0=[0.8 ; 1.2];

v=[v0]; V=[v0];



for i=1:1:15
    
    T=t(i);
    
    funcioP5=@(x)([log((3*x(2)-1)/(3*x(1)-1))+9/(4*T)*(1/x(2)-1/x(1))-1/(3*x(2)-1)+1/(3*x(1)-1) (8*T)/3*(1/(3*x(2)-1)-1/(3*x(1)-1))-1/x(2)^2+1/x(1)^2]);
    
    [XK,resd,it] = newtonn(v,10^-16,100,funcioP5);
    
    v=XK(end);
    
    V=[V' , XK(end)];
    
    
end
