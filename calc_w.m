function [W,U] = calc_w(p,x,n,y);

% -- [W,U] = calc_w(p,x,n);
% -- [W,U] = calc_w(p,x,n,y);
%
%
% The purpose of this function is to calculate both the W
% matrix and U matrix in a given system, either for the
% resident only (no y input) or for the mutant (specify
% mutant traits y as well). 
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
% n: Population size. Row vector.
%
% y: The hatching-date strategy being pursued by the
% variant ("mutant strategy"). Row vector. Optional.
%
%
% OUTPUTS
%
% W: The recruitment rate matrix.
%
% U: The growth rate matrix.
%


h = length(x);

if nargin == 3
    % Doing resident
    E = calc_E(p,x);
    F = calc_F(p,x,n);
else
    % Doing mutant
    E = calc_E(p,y);
    F = calc_F(p,x,n,y);
end

A = max(0,(p.K-p.s.*n)./p.K);

pM = p.pM; % Construct pM matrix of dispersal

% Construct W matrix: W(i,j) is the flow from hab i to hab j
W = repmat(E',1,h).*repmat(p.s,h,1).*F.*repmat(A,h,1);
W = pM.*W; % Put this here to match Jacob's for debugging

if length(p.s) == 1
    U = (W+p.s*eye(h,h))'; % Note transposition
else
    U = (W+diag(p.s))';
end
