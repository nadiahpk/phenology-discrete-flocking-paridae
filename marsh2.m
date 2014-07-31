% For Marsh Tits as described in Nilsson & Smith (1988) , but with habitats set equal. 

p.K = [50 50]; 
p.m = 0.8; 
p.disptype = 2; % 2 = Tony's way

p.x_opt = [131,171]; % It's the relative distance between these that matters

p.E_0 = [7 7]; % Wikipedia, clutch of 5-9 eggs, modified to make equal
p.sigma_E = 1.5*[6 6]; % See supplementary materials to paper

% v = 193.62 - from paper
% m = 15; % Choose intermediate value for illustrative purposes
p.mu_g =  2.3976;
p.sigma_g =  0.78795;

% Mean and standard dev of underlying normal distribution
p.v_g = exp(2*p.mu_g+p.sigma_g^2)*(exp(p.sigma_g^2)-1);
p.m_g = exp(p.mu_g+p.sigma_g^2/2);

p.s = [0.6 0.6]; % Arbitrary, usually about a half

% Integration parameters
p.Nfn = 1; % Number of grid points for calculating f. Can set to 1.
p.Ncrs = 100; % Number of grid points for calculating F
p.tol = 0.0001; % End integration after the cumulative arrival distribution of last to arrive less than this value
p.adaptsimp = 1e-5; % Tolerance on adaptive simpson's method to determine integration steps. 1e-5 is adequate for most except when we get to quite small mg (mg ~< 7).

% Differentiation parameters
p.del = 1e-4;

% Two ways to calculate pM, type "2" used in paper 
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

