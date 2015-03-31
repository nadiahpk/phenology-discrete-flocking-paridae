function [eigH,eigJ,eigJs,Hess,Jac] = check_stab(p,x,n)

% -- [eigH,eigJ,eigJs,Hess,Jac] = check_stab(p,x,n)
%
% The purpose of this function is to check the evolutionary
% stability of the singular strategy. This is done by calculating the
% eigenvalues of the Hessian matrix and the symmetric
% Jacobian of the fitness gradient. If the eigenvalue is
% negative the singular strategy passes the corresponding
% test (below).
%
% INPUTS
%
% p: Dictionary of parameter values
%
% x: Trait values at evolutionarily singular strategy
%
% n: Population sizes at evolutionarily singular strategy
%
% OUTPUTS
% 
% eigH: Dominant eigenvalue of the Hessian
% 
% eigJ: Dominant eigenvalue of the Jacobian
% 
% eigJs: Dominant eigenvalue of the J^s (Leimar 2009)
% 
% Hess: The Hessian matrix
% 
% Jac: The Jacobian matrix


% Modify integration and derivative parameters for numerical
% solutions:
% 1. Improve the integration 
p.adaptsimp = p.adaptsimp/10; 
% 2. Double derivative step sizes need to be lower than that
% used to find the singular strategy
p.del = 5000*p.del; % == 0.5 for default
del = 2*p.del; % == 1

% Prepare storage
Hess = [];
Jac = [];

% Calculate second derivatives
for ind = 1:length(x) % For each trait

    yhi = x; ylo = x;

    yhi(ind) = yhi(ind)+del/2;
    ylo(ind) = ylo(ind)-del/2;

    g_yhi = calc_dfit(p,yhi',n,x',1); % no dens dep
    g_ylo = calc_dfit(p,ylo',n,x',1);
    g_dy = (g_yhi-g_ylo)/del;

    Hess = [Hess,g_dy];

    g_yhi = calc_dfit(p,yhi',n,yhi',0); % dens dep
    g_ylo = calc_dfit(p,ylo',n,ylo',0);
    g_dy = (g_yhi-g_ylo)/del;

    Jac = [Jac,g_dy];

end

% Eigenvalues for checking stability
eigH = max(real(eig(Hess)));
eigJ = max(real(eig(Jac)));
Js = (Jac+Jac')/2;
eigJs = max(real(eig(Js)));
