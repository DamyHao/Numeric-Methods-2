function UR = qrFact(A)
    [n, m] = size(A);
    if()
   
%{
  if (m > n)
    
    end 
%}

M = n;
for i = 1:1:M
    B = A(i:M, i:M); % Es la A petita
    x = B(:, 1); % Es la primera columna de la A petita
    tau = sign(x(1))*norm(x); 
    gama = 1 + x(1)/tau;
    u = [1; x(2:end)/(x(1)+tau)];
    vAux = gama*u'*B;
    A(i:M, i:M) = B - u*vAux; 
end
