function neqm = calc_n(p,x,n0);

% -- neqm = calc_n(p,x,n0);
% -- neqm = calc_n(p,x);
%
%
% The purpose of this function is to calculate the
% steady-state population size given the parameters and
% hatching-date strategy.
%
%
% INPUTS
%
% p: The dictionary of parameter values. See marsh2.m for an
% example of how to specify these.
%
% x: The hatching-date strategy being pursued by the
% resident species. Row vector.
%
% n0: An initial estimate for the steady-state population size
% at the evolutionarily singular strategy. Row vector.
%
%
% OUTPUTS
%
% neqm: The steady state population size. Row vector.
% Optional.
%

if nargin == 2
    n0 = p.K; % Need to find a better way of doing this
end

f = @(n) calc_deln(p,x,n);
[neqm,err] = fsolve(f,n0');
neqm = neqm';
