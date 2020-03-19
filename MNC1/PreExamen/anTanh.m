function I=tanH(a,b,n,f,C)

h=C/(n)^.5;

sumatori=[];
for j=-n:n
   
    xj=(b-a)*.5*tanh(j.*h/2)+(b+a)*.5;
    fj=f(xj);
    sumatori=[sumatori fj.*(cosh(j.*h/2)).^-2];
    
end

s=sum(sumatori);
I=.5*h*s*(b-a)*.5;
end