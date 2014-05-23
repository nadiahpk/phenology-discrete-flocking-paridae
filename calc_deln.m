function deln = calc_deln(p,x,n0);

% Accepts a vertical n0 from solver, so some confusing transposing is
% necessary

[W,U] = calc_w(p,x,n0');
n1 = U*n0;
%n1 = min([n1';p.K])'; % Shouldn't need to do this
deln = n1 - n0;


