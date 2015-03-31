function p = updated_pm(p);

% -- p = upated_pm(p)
%
%
% The purpose of this function is to update both the arrival
% distribution parameters (p.mu_g and p.sigma_g) and the
% dispersal matrix (p.pM) in the parameters dictionary p.
%
% In many cases we would like to explore the effect of
% changing a parameter value on various behaviours, however
% there is interdependence between parameters within the
% parameters dictionary p. For example, if p.m_g changes
% then p.mu_g and p.sigma_g will change as well. Therefore
% these will need to be updated.
%
%
% INPUTS 
%
% p: The dictionary of parameter values. See
% default_parameter_values.m for an example of how to
% specify these.
%
%
% OUTPUTS
%
% p: A consistent dictionary of parameter values. 

% Two ways to calculate pM 
if p.disptype == 1; % Jacob's way
    h = length(p.x_opt);
    if h == 1;
        % If only one habitat stay where you are
        p.pM = 1; 
    else
        % Disperse with probability m to h-1 other habitats
        % equally
        p.pM = (p.m/(h-1))*(ones(h,h)-eye(h,h)) + (1-p.m)*eye(h,h);
    end
end
if p.disptype == 2; % Tony's way
    h = length(p.x_opt);
    if h == 1;
        % If only one habitat stay where you are
        p.pM = 1; 
    else
        % Disperse with probability m to h-1 other habitats
        % equally
        p.pM = (1-p.m)*eye(h,h) + p.m*repmat(p.K/sum(p.K),h,1);
    end
end

% Changes to arrival distribution
p.mu_g = log(p.m_g^2/(sqrt(p.v_g+p.m_g^2)));
p.sigma_g = sqrt(log(1+p.v_g/p.m_g^2));
