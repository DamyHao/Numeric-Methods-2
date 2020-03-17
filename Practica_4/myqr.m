function UR = myqr(A)
    % qr del alvaro
    [n, m] = size(A);
    if n > m; M = m; else M = m - 1; end

    for jj = 1:M
        B = A(jj:n, jj:m);
        x = B(:, 1); tau = sign(x(1)) * norm(x); gam = 1 + x(1) / tau;
        u = [1; x(2:end) / (tau + x(1))]; vauxT = gam * u' * B;
        A(jj:n, jj:m) = B - u * vauxT; A(jj + 1:n, jj) = u(2:end);
    end

    UR = A;
