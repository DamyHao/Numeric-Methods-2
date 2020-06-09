function DX=fun2(X,t)
x=X(1,1);px=X(2,1);
y=X(3,1);py=X(4,1);
dx  = px;
dpx = -4*(x^2+y-11)*x -2*(x+y^2-7);
dy  = py;
dpy = -2*(x^2+y-11) -4*(x+y^2-7)*y;
DX=[dx;dpx;dy;dpy];
end