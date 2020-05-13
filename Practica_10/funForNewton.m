function r = funForNewton(z)
    % Function to launch newtonn
    % Input: z is a vector of two components, the first one corresponds to
    % the launch angle of the particle and the second one to the time that it
    % takes to get to the origin
    % Output: the distance to the origin. Newton will make it 0.
    initial = [0; 0; sqrt(2) .* cos(z(1)); sqrt(2) .* sin(z(1))]; % initial point to launch RK4
    steps = 20000; %h will be smaller than 0.0001 since we the time to be < 2

    if 0 < z(2) < 2
        h = z(2) / steps;
        sol = RK4(initial, h, @gravFunctionV, steps + 1);
        r = sol(1:2, end);
    else % If time is greater than 2, we will introduce an "artificial slope" to help newton to converge to r = (0,0)
         % If we launch newton to the correct point we will not reach this code:
        disp('Surpasing t = 2');
        r = [1; 1]*z(2);
    end
end