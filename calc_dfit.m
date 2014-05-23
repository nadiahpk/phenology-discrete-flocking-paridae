function df = calc_dfit(p,y,n0,x);

% INPUTS
%
% y: Variant strategy
%
% x: Resident strategy

y = y'; % Accepts y as a column vector, for use with fsolve

if nargin < 4
    x = y;
    if nargin < 3
        n0 = p.K;
    end
else
    x = x'; % Accepts x as a column vector, for use with fsolve
end
n = calc_n(p,x,n0); % Calculates ss for resident trait

del = p.del;
df = [];

for i = 1:length(y) % For each trait/habitat type

    yd = y;
    yd(i) += del/2;
    fit_hi = calc_fit(p,x,n,yd);

    yd = y;
    yd(i) -= del/2;
    fit_lo = calc_fit(p,x,n,yd);

    df = [df; (fit_hi-fit_lo)/del];

end
