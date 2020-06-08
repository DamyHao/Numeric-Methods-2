% BDF1 (implicit Euler time-stepper)
% Solution for u_t = f(t,u) with u(t0) = v0 (n-dimensional)
% INPUT
%   fun: funcio a integrar. Ha de ser del tipus fBurger=@(t,f)(0.1*D2*f+f.*D1*f); (encara que no s'utilitzi t)
%   t0: (initial time)
%   h: (time-step) ; v0 (initial condition) ; N (no. steps)
% OUTPUT:
%   T: (time vector: 1 x N+1)
%   Y: (solution matrix: n x N+1)
function [T,Y] = bdf1(fun,t0,h,v0,N)
T = zeros(1,N+1); n = length(v0); I = eye(n);
Y = zeros(length(v0),N+1); 
DF = zeros(n); % El que sera el Jacob
H = sqrt(eps)*I; % Matriu per fer les derivades direccionals del Jacobia
T(1) = t0; Y(:,1) = v0;
for j = 1:N
    z0 = Y(:,j) ; tplus = t0 + h*j; T(j+1) = tplus;
    dz = 1; z = z0 ;
    % Aqui es com si comences newton
    while norm(dz) > 1e-12
        f1 = feval(fun,tplus,z);
        % Fa el Jacobia:
        for kk = 1:n
            f2 = feval(fun,tplus,z+H(:,kk));
            Df(:,kk) = (f2 - f1)/H(kk,kk);
        end
        % Equacio que estem intentant resoldre (per aixo es implicit)
        %v_n+1 = z, vn = x0 f(t_n+1, v_n+1) = f1
        % Volem resoldre per quion valor de z es compleix z-z0-h*f1 == 0
        dz = -(I-h*Df)\(z-z0-h*f1); z = z + dz;
    end
    Y(:,j+1) = z ; % Trobada ya la solucio amb newton, podem establir el seguente element
    % de la successio.
end