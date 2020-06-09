function [dX]=funb(t,X)
x =X(1,1);
px=X(2,1);
y =X(3,1);
py=X(4,1);

dx=px;
dpx=-4*x*(x^2+y-11)-2*(x+y^2-7);
dy=py;
dpy=-4*y*(x+y^2-7)-2*(x^2+y-11);
dX=[dx;dpx;dy;dpy];
end