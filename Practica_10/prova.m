n = -3:0.01:-1;
haches = 10.^n
%{
 
for h = haches
    points = time / h + 1;
    %RK i AB gives us the points specified in the arguments so to obatin the one
    %corresponding to t=2 we have to add this plus one to the steps

    if floor(points) == points
        % The number of points must be a natural number and for some values
        %of h it is a decimal one, we one use those values that gives a
        %natural number of points
        hs = [hs h];
        solutionRK = RK4(initial, h, @gravFunctionV, points);
        r = solutionRK(1:2, end);
        errorsRK = [errorsRK norm(r - rext)];

        %As done before, we use the first three steps provided by RK plus
        %the intial condition
        solutionAB = AB4(solutionRK(:, 1:4), h, @gravFunctionV, points);
        r2 = solutionAB(1:2, end);
        errorsAB = [errorsAB norm(r2 - rext)];
    end
end 
%}