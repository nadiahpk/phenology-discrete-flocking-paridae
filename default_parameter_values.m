% Default parameter values, see "Appendix S1: Default model
% parameter values" for a discussion of how they were
% selected

% Tip: These are good start values for this parameter-value set
% x0 = [162.91   162.90]; n0 = [46.379   46.378];


% Parameter values from the literature
% -------------------------------------

% Survival probability, based on that reported for willow
% tits in Lahti et al. 1996
p.s = [0.6 0.6]; 

% Arrival distribution.
% For Marsh Tits as described in Nilsson & Smith (1988), I
% estimate p.mu_g = 1.6 and p.sigma_g = 1.1. This gives m_g
% = 9.073 and v_g = 193.62. Marsh tits are probably near low
% end of mean arrival time, so we made the default m_g
% higher. Keep v = 193.62 and make m_g = 15. Therefore
p.mu_g =  2.3976;
p.sigma_g =  0.78795;


% Maximum reproduction rate. A clutch size of 5-9 eggs is
% usual. Lambrechts et al. (1997) gives values of
% approximately 7.2 and 4.6 fledglings per nest in Corsican
% deciduous and evergreen habitats respectively. I set the
% default at the high end of reproduction, in order to show
% the phenological effect more clearly.
p.E_0 = [4.5 4.5]; 


% Width of the relationship between food availability and
% hatching time. No direct information about this exists in the
% literature, however Vatka et al. (2011) reports an
% asynchrony of approximately 10 days. Therefore, find
% p.sigma_E such that when an arrival distribution with mu_g
% = 1.6, sigma_g = 1.1 is used, it results in approximately
% 10-day asynchrony
p.sigma_E = [10 10]; 


% Parameter values that depend upon the scenario, explore range
% -------------------------------------------------------------

% Landscape composition [early habitat %, late habitat %]
p.K = [50 50]; 

% Dispersivity (Equation 2)
p.m = 0.8; 
p.disptype = 2; % Dispersal model style, 2 = Tony Ives' way of modelling dispersal.

% Optimal hatching time in each habitat, [\hat{x} in early habitat, \hat{x} in late habitat] 
% It's the difference between these that matters in this version of the model
p.x_opt = [171,171]; 



% Integration and differentiation parameters
% ------------------------------------------
% One shouldn't need to change these

p.Nfn = 1; % Number of grid points for calculating f. Can set to 1.
p.Ncrs = 100; % Number of grid points for calculating F
p.tol = 0.0001; % End integration after the cumulative arrival distribution of last to arrive less than this value
p.adaptsimp = 1e-5; % Tolerance on adaptive simpson's method to determine integration steps. 1e-5 is adequate for most except when we get to quite small mg (mg ~< 7).

% Differentiation parameters
p.del = 1e-4;


% Parameter values derived from other parameter values
% ----------------------------------------------------

% Mean and standard dev of underlying non-logarithmised
p.v_g = exp(2*p.mu_g+p.sigma_g^2)*(exp(p.sigma_g^2)-1);
p.m_g = exp(p.mu_g+p.sigma_g^2/2);

% Dispersal: Two ways to calculate pM. Type "2" is used in
% paper.
if p.disptype == 1; 
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
if p.disptype == 2; 
    h = length(p.x_opt);
    if h == 1;
        % If only one habitat stay where you are
        p.pM = 1; 
    else
        % Disperse proportional to availability
        p.pM = (1-p.m)*eye(h,h) + p.m*repmat(p.K/sum(p.K),h,1);
    end
end


