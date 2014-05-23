function [W,U] = calc_w(p,x,n,y);

% Calculates W matrix either for resident (no y input) or
% for mutant (specify mutant traits y)
% NOTE: You'll need to transpose the W before using it in
% pop dynamic equation

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

P = max(0,(p.K-p.s.*n)./p.K);

pM = p.pM; % Construct pM matrix of dispersal

% Construct W matrix: W(i,j) is the flow from hab i to hab j
W = repmat(E',1,h).*repmat(p.s,h,1).*F.*repmat(P,h,1);
W = pM.*W; % Put this here to match Jacob's for debugging

if length(p.s) == 1
    U = (W+p.s*eye(h,h))'; % Note transposition
else
    U = (W+diag(p.s))';
end
