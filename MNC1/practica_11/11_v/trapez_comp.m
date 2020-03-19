function I_trap_c = trapez_comp(a, b, m, fun)

k = 0:m;
H = (b-a)/m;
xk = a+k*H;
sumatori = 0;
for i=0:m-1
    sumatori = sumatori + (fun(xk(i+1))+fun(xk(i+2)));
end
I_trap_c = (H/2)*sumatori;
end

