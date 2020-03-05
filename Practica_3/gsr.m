function [Q,R] = gsr(A)

[n, m] = size(A);
R = zeros(m);
Q = 0 * A;

R(1, 1) = norm(A(:, 1)); Q(:, 1) = A(:, 1) / R(1, 1);
S = zeros(m);

for k = 2:m
    R(1:k - 1, k) = Q(:, 1:k - 1)' * A(:, k); %No es k-1

    
    Q(:, k) = A(:, k) - Q(:, 1:k - 1) * R(1:k - 1, k);
    
    S(1:k-1,k) = Q(:,1:k-1)'*Q(:, k); % matriu que conté els productes <q_j,q(tilde)_k)
    R(k, k) = norm(Q(:, k)); % Definim Rkk com a teoria.
    
   
    Q(:, k) = Q(:, k) / R(k, k); % Normalitzem el vector q.
end


end
