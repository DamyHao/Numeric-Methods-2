% Function that performs Newton method to solve n-dimensional nonlinear systems of equations
% Input:  x0 - [Rn column vector] - Point in Rn, initial guess of the solution to launch the method, must be
%              close to the desired solution.
%        tol - [scalar] - Desired tolerance, once the solutions differ less
%              than the tolerance, the method stops.
%      itmax - [scalar] - Maximum number of iterations, the method stops
%              after itmax iterations.
%          f - [vector function (Rn --> Rn)] - Function that gives rise to
%              the system f(x)=0 that we want to solve.
%
% Output:  XK - [n x it matrix] - kth column of XK contains the kth aproximation to
%               the solution
%         res - [Rn column vector] - kth component of res contains the
%               residual ||f(xk)|| (Should be 0 in the exact solution).
%          it - [scalar] - Number of iterations performed.
% Note: Function modified to use a parameter T in the evaluations of f
% Requires jacmod.m, factPPLU.m and solvePPLU.m

function [XK,res,it]=newtonNLSmod(x0,tol,itmax,f,T) 
xk=[x0]; res=[norm(f(xk,T))]; XK=[xk];it=1;
tolk=1; 
while it < itmax & tolk > tol
    fk=f(xk,T); Jf=jacmod(f,xk,T);
    [P,L,U]=factPPLU(Jf);
    dx=solvePPLU(L,U,P,-fk);
    xk=xk+dx;
    XK=[XK xk]; res=[res norm(fk)];
    tolk=norm(XK(:,end)-XK(:,end-1));
    it=it+1;
end