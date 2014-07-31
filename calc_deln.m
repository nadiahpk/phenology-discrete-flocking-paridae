function deln = calc_deln(p,x,n0);

% -- deln = calc_deln(p,x,n0);
%
%
% The purpose of this function is to calculate the rate of
% change of the population size. It is used by calc_n.m to
% find the population steady-state. Note that because it is
% used with fsolve it accepts the n0 as a column vector.
%
%
% INPUTS
%
% p: The dictionary of parameter values. See marsh2.m for an
% example of how to specify these.
%
% x: The hatching-date strategy being pursued by the
% population ("resident strategy"). Row vector.
%
% n0: The population size at time 0. Column vector.
%
%
% OUTPUTS
%
% deln: The derivative of the population size.


[W,U] = calc_w(p,x,n0');
n1 = U*n0;
%n1 = min([n1';p.K])'; % Shouldn't need to do this
deln = n1 - n0;


