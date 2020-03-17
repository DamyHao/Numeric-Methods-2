% Code 18: GMRES (Arnoldi iteration - Reorthogonalized)
function [x, k] = mygmres(Afun, b, tol, dimkryl)
    % Requires: myqr.m and qrsolve.m (Codes 15 & 16)
    % Input: Afun (x --> Ax function) | b (r.h.s)
    % tol (iteration stops if |A*x_m-b| < tol)
    % dimkryl: max. dimension of krylov space.
    % Output: x: approximation found
    % k: no. of iterations used
    hkp1 = 1; resd = 1; k = 1;
    Q(:, 1) = b / norm(b); H = []; % q_1 = b/|b|

    while resd > tol & hkp1 > eps & k < dimkryl
        Q(:, k + 1) = feval(Afun, Q(:, k)); % q_k+1 = Aq_k
        h = Q(:, 1:k)' * Q(:, k + 1); % h_{ij} - coeffs.
        Q(:, k + 1) = Q(:, k + 1) - Q(:, 1:k) * h; % See (5.255)
        S = Q(:, 1:k)' * Q(:, k + 1); Q(:, k + 1) = Q(:, k + 1) - Q(:, 1:k) * S;
        h = h + S; % Reorth.
        hkp1 = norm(Q(:, k + 1)); H(1:k + 1, k) = [h; hkp1];
        Q(:, k + 1) = Q(:, k + 1) / hkp1;
        UR = myqr(H(1:k + 1, 1:k)); e1b = [norm(b); zeros(k, 1)];
        y = qrsolve(UR, e1b); x = Q(:, 1:k) * y'; k = k + 1;
        resd = norm(feval(Afun, x) - b);
    end
