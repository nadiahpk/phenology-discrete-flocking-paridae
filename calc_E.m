function E = calc_E(p,x)

% -- E = calc_E(p,x)
%
%
% The purpose of this function is to calculate the
% reproductive output in each habitat type given the
% hatching-date strategies pursued.
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
%
% OUTPUTS
%
% E: Reproduction rate. Row vector.

E = p.E_0.*exp(-(x-p.x_opt).^2./(2*p.sigma_E.^2));

