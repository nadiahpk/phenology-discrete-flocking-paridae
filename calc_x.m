function x = calc_x(p,x0,n0);

% Calculate singular hatching time strategy

if nargin == 1;
    % Initialise with local optima and max poss. populations
    x0 = p.x_opt;
    n0 = p.K;
end

% Steady state
f = @(u) calc_dfit(p,u,n0);
[x,errs] = fsolve(f,x0');
x = x'; % return as row vector
