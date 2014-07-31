function x = calc_x(p,x0,n0);

% -- x = calc_x(p,x0,n0);
%
%
% The purpose of this function is to calculate the
% evolutionarily singular hatching-time strategy. This
% function will make an attempt if given no starting-point
% values, but from experience it's a good idea to give it a
% reasonably good first estimate if you want it to find the
% singular point in reasonable time (or at all).
%
%
% INPUTS
%
% p: The dictionary of parameter values. See marsh2.m for an
% example of how to specify these.
%
% x0: An initial estimate for the singular strategy. Row
% vector.
%
% n0: An initial estimate for the steady-state population size
% at the evolutionarily singular strategy. Row vector.
%
%
% OUTPUTS
%
% x: A numerical estimate for the evolutionarily singular
% strategy. Row vector.
%

if nargin == 1;
    % Initialise with local optima and max poss. populations
    x0 = p.x_opt;
    n0 = p.K;
end


% Steady state
f = @(u) calc_dfit(p,u,n0);
[x,errs] = fsolve(f,x0');
x = x'; % return as row vector
