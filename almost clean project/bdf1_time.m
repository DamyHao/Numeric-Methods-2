% Code 27: BDF1 (implicit Euler time-stepper)
% Solution for u_t = f(t,u) with u(t0) = v0 (n-dimensional)
% Input: fun (function name) ; t0 (initial time)
% h (time-step) ; v0 (initial condition) ; N (no. of points)
% Output: T (time vector: 1 x N+1)
% Y (solution matrix: n x N+1)


%AQUEST BDF1 NECESSITA QUE LA FUN QUE LI ENTREM TINGUI UNA VARIABLE
%D'ENTRADA QUE SIGUI EL TEMPS, SI NO ES EL CAS UTILITZEM BDF1.
%LA FUNCIÓ D'ENTRADA, EN EL PRIMER INPUT HA DE TENIR LA VARIABLE TEMPS, QUE
%SI NO ES UTLITIZADA S'HA DEXPRESAR AMB ~, EXEMPLE:

%{
function v=funcioexemple(~,x)
v(1)=x(2);
v(2)=50*(1-x(1)^2)*x(2)-x(1);

v=[v(1);v(2)];
end
%}

function [T,Y] = bdf1_time(fun,t0,h,v0,N)
T = zeros(1,N); n = length(v0); I = eye(n);
Y = zeros(length(v0),N); DF = zeros(n); H = sqrt(eps)*I;
T(1) = t0; Y(:,1) = v0;
for j = 1:N-1
    z0 = Y(:,j) ; tplus = t0 + h*j; T(j+1) = tplus;
    dz = 1; z = z0 ;
    while norm(dz) > 1e-12
        f1 = feval(fun,tplus,z);
        for kk = 1:n
            f2 = feval(fun,tplus,z+H(:,kk));
            Df(:,kk) = (f2 - f1)/H(kk,kk);
        end
        dz = -(I-h*Df)\(z-z0-h*f1); z = z + dz;
    end
    Y(:,j+1) = z ;
end