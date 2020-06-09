%f: must be f(t,X) with X column vector;
%N: number of time steps
%h: time distance between steps
%v0:initial point

function [v,t]=RK4(v0,N,h,f)
v=v0; t=(0:h:h*N);
for i=2:length(t)
    A=h*f(t(i),v(:,end));
    B=h*f(t(i)+h/2,v(:,end)+A/2);
    C=h*f(t(i)+h/2,v(:,end)+B/2);
    D=h*f(t(i)+h,v(:,end)+C);
    vn=v(:,end)+(1/6)*(A+2*B+2*C+D);
    v=[v vn]; 
end

