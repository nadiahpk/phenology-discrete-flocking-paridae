function [f,gg] = calc_fg(p,x,n,x_crs);

% -- [f,gg] = calc_fg(p,x,n,x_crs);
%
%
% The purpose of this function is to return the probability
% of arriving from each habitat (gg) at each date corresponding
% to the dates specified in x_crs, and the probability of
% being first in the queue (f), determined by integrating on
% the grid specified by x_crs, at each date x_crs. Its
% primary use is in calc_F, which is where the overall
% probability of being first in the queue is calculated. The
% integration is performed using the trapezoid method.
%
%
% INPUTS
%
% p: The dictionary of parameter values. See marsh2.m for an
% example of how to specify these.
%
% x: The hatching-date strategy being pursued by the
% population. Row vector.
%
% n: The population size. Row vector.
%
% x_crs: A grid of arrival dates upon which the arrival
% distribution is evaluated and the integration to evaluate
% f is performed.
%
%
% OUTPUTS
%
% f: The probabilities of being first in a territory queue
% given the arrival dates in the corresponding x_crs.
%
% gg: The probabilities of arrival on the dates
% corresponding x_crs.
%

h = length(x); % Number of habitats, traits
pM = p.pM; % Dispersal probability matrix

Nfn = p.Nfn;
Ncrs = length(x_crs);

if Nfn == 1;
    % No integral to calculate f at each course grid point
    % Surprisingly, this still gives fairly accurate results, and it's
    % a good deal quicker
    del = x_crs(2)-x_crs(1);
    gg = lognpdf(repmat(x_crs',1,h)-repmat(x,Ncrs,1),p.mu_g,p.sigma_g);
    r = repmat(n.*calc_E(p,x),Ncrs,1).*gg;
    lam = r*pM.*repmat(p.s,Ncrs,1)./repmat(p.K,Ncrs,1);
    %cum_lam = del*cumtrapz(lam);
    cum_lam = cumtrapz(x_crs',lam); % An alternative for unevenly spaced 
else
    % Do first point
    cum_lam = zeros(1,h); % Integral from x_lo to x_lo = 0
    gg = lognpdf(x_lo*ones(1,h)-x,p.mu_g,p.sigma_g);

    % Do remainder
    for crs_ind = 2:Ncrs
        x_crs_lo = x_crs(crs_ind-1);
        x_crs_hi = x_crs(crs_ind);
        x_fn = linspace(x_crs_lo,x_crs_hi,Nfn); % Grid within course grid

        % Calculate partial integral of lambda
        g = lognpdf(repmat(x_fn',1,h)-repmat(x,Nfn,1),p.mu_g,p.sigma_g);
        r = repmat(n.*calc_E(p,x),Nfn,1).*g;
        lam = r*pM.*repmat(p.s,Nfn,1)./repmat(p.K,Nfn,1);
        part_lam = (x_fn(2)-x_fn(1))*trapz(lam);

        % Add to cum_lam vector
        cum_lam(crs_ind,:) = cum_lam(crs_ind-1,:) + part_lam;

        % Store the gg while you're at it
        gg(crs_ind,:) = lognpdf(x_crs_hi*ones(1,h)-x,p.mu_g,p.sigma_g);
    end
end

% Calculate f on course grid
f = exp(-cum_lam);

