y0 = [-4 ; 4] ; y1 = [-3.975 ;4.03282550110059 ];
xv = [y0(1) y1(1)] ; yv = [y0(2) y1(2)] ; s = 1 ;
for ii = 1:1000
[y,iconv] = continuationStep(@fun,y0,y1,s,1e-12,10);
xv = [xv y(1)] ; yv = [yv y(2)] ; y0 = y1 ; y1 = y ;
end
plot(xv,yv,'-k'); axis equal ; axis([-5.5 5.5 -5.5 5.5])