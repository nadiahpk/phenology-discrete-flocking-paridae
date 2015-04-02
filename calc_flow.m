function flow = calc_flow(p,x,n);

% -- flow = calc_flow(p,x,n);
% -- flow = calc_flow(p);
%
%
% Quantify the flow of new recruits into each habitat from
% which habitat
%
%
% INPUTS
%
% p: The dictionary of parameter values. See
% default_parameter_values.m for an
% example of how to specify these.
%
% x: The hatching-date strategy.
%
% n: The population size at that hatching-date strategy.
%
%
% OUTPUT
%
% flow: A matrix showing the proportion of recruits in each
% habitat that come from each habitat. flow(i,j) is the
% proportion of recruits in habitat i that came from habitat
% j.

if nargin < 2
    x = calc_x(p);
    n = calc_n(p,x);
    x = calc_x(p,x,n);
    n = calc_n(p,x,n);
end

[W,U] = calc_w(p,x,n);
h = length(n);
flow = W.*repmat(n',1,h);
flow = flow./repmat(sum(flow,1),h,1); % Normalises it so habitats total 1
