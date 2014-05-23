function neqm = calc_n(p,x,n0);

% Accepts a horizontal n0

if nargin == 2
    n0 = p.K; % Need to find a better way of doing this
end

f = @(n) calc_deln(p,x,n);
[neqm,err] = fsolve(f,n0');
neqm = neqm';
