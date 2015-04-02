function fitness = calc_fit(p,x,n,y);

% -- fitness = calc_fit(p,x,n,y);
%
%
% The purpose of this function is to calculate the invasion
% fitness of a variant strategy y in a population following
% the resident strategy x. 
%
%
% INPUTS
%
% p: The dictionary of parameter values. See
% default_parameter_values.m for an
% example of how to specify these.
%
% x: The hatching-date strategy being pursued by the
% population ("resident strategy"). Row vector.
%
% n: The steady-state population size. Row vector.
%
% y: The hatching-date strategy being pursued by the
% variant ("mutant strategy"). Row vector. Optional.
%
%
% OUTPUTS
%
% fitness: Invasion fitness of the variant strategy. A value
% greater than 1 means that the variant can invade.
%

h = length(x);

if nargin == 3
    [W,U] = calc_w(p,x,n);
else
    [W,U] = calc_w(p,x,n,y);
end

fitness = max(real(eig(U)));
