function df = calc_dfit(p,y,nn0,x,n_flag);

% -- df = calc_dfit(p,y,nn0,x,n_flag);
% -- df = calc_dfit(p,y,nn0);
%
%
% The purpose of this function is the calculate the
% fitness gradient. Note that, because of it is designed to
% be used with fsolve, it accepts y (and x) as column vectors.
%
% This function is used in two ways. In calc_x.m
% it is used to find the singular strategy by evaluating the
% slope around the resident strategy (n_flag = 0), and in
% check_stab.m it is used to approximate the second
% derivative of the slope by evaluating it at different mutant
% strategies.
%
%
% INPUTS
%
% p: The dictionary of parameter values. See marsh2.m for an
% example of how to specify these.
%
% y: The hatching-date strategy being pursued by the
% variant ("mutant strategy"). Column vector. Optional.
%
% nn0: Either the initial guess for the steady-state
% population size (n_flag = 0) or the actual steady-state
% population size (n_flag = 1). Row vector.
%
% x: The hatching-date strategy being pursued by the
% population ("resident strategy"). Column vector.
%
% n_flag: When set to 1, indicates that nn0 represents the
% steady-state population size, otherwise it is treated as
% an initial guess for the steady-state and the steady-state
% is solved with calc_n.m.
%
%
% OUTPUTS
%
% df: The fitness gradient. When this is zero there is an
% evolutionarily singular strategy. Column vector.

y = y'; % Accepts y as a column vector, for use with fsolve

if nargin < 5
    n_flag = 0;
end

if nargin < 4
    x = y;
    if nargin < 3
        nn0 = p.K;
    end
else
    x = x'; % Accepts x as a column vector, for use with fsolve
end

if n_flag == 1
    n = nn0;
else
    n = calc_n(p,x,nn0); % Calculates ss for resident trait
end

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
