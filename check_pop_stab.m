function eigPopJ = check_pop_stab(p,x,n)

% -- eigPopJ = check_stab(p,x,n)
%
% The purpose of this function is to check the population-dynamic
% stability of the singular strategy. This is done by linearising
% the system and then numerically checking that the dominant 
% eigenvalue of the Jacobian has absolute value less than 1.
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
% eigPopJ: Dominant eigenvalue of the Jacobian of the pop
% dynamics

h = length(x);
del = 1; % perturbation size we're using
J = zeros(h,h); % place to store the Jacobian

for j = 1:h

    pert = n; pert(1,j) = pert(1,j) + del;
    [W,U] = calc_w(p, x, pert);

    % numerical estimate of derivative
    dfi_dfj = (U*pert' - n')/del;

    % append to jacobian
    J(:,j) = dfi_dfj;

end

eigPopJ = max(eig(J));
