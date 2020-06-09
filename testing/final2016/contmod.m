function [y,iconv]=contmod(f,y0,y1,s,tol,itmax)
tolk=1; it=0; n=length(y0)-1;
v=y1-y0; yp=y1+s*v; xk=yp;
while it < itmax & tolk > tol
    Fk=[f(xk); v'*(xk-yp)]; DFk=[jac(f,xk); v'];
    [P,L,U]=factPPLU(DFk); dxk=solvePPLU(L,U,P,-Fk);
    xk=xk+dxk; tolk=norm(dxk); it=it+1;
end
y=xk;
if it<=itmax & tolk<tol
    iconv=0;
else
    iconv=1;
end