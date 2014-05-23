% For Marsh Tits as described in Nilsson & Smith (1988) , but with habitats set equal. 

p.K = [50 50]; % Equal
p.m = 1; % guess
p.disptype = 2; % 2 = Tony's way

%p.htod = 32; % Number of days from hatching to natal dispersal
%p.x_opt = [141,141+30]; % x1 is deciduous from Nilsson, and x2 is approximately when evergreen peak is
p.x_opt = [131,171]; % Both set to the evergreen peak from Nilsson

p.E_0 = [7 7]; % Wikipedia, clutch of 5-9 eggs, modified to make equal
p.sigma_E = 1.5*[6 6]; % A bit of a guess, based off Veen et al. (2010), I selected the 1.5 so that, when the code is run using the arrival distribution below, the hatching date selected is about 10-15 days before the optimal hatching date. I got the 10-15 days from Vatka et al. (2011) and Both et al. (2008).

% See 140404 for how these parameters were found. The short of it is that they were selected for agreement with Nilsson and Smith (1988)
% m = 9.073, v = 193.62 - from paper
%p.mu_g = 1.6; 
%p.sigma_g = 1.1;
% The solution to this system is near:
% x0 = [159.25   159.25]; n0 = [45.674   45.674];
% Note that the solver has difficulty finding this point without a good starting point

% Choose a slightly less painful intermediate value to use
% as the default for exploring the parameter space. This
% value also gives us more variety in the parameter
% response.
% m = 15; v as above
p.mu_g =  2.3976;
p.sigma_g =  0.78795;
% The solution to this system is near:
%   x0 = [124.96   164.85]; n0 = [4.9258e+01      4.9258e+01];
% or for p.x_opt = [171 171]:
%   x0 = [161.88   161.88]; n0 = [48.155   48.155];

% Arbitrary, usually about a half
p.s = [0.6 0.6]; 

% Integration parameters
p.Nfn = 1; % Number of grid points for calculating f. Can set to 1.
p.Ncrs = 100; % Number of grid points for calculating F
p.tol = 0.0001; % End integration after the cumulative arrival distribution of last to arrive less than this value

% Differentiation parameters
p.del = 1e-4;

% Two ways to calculate pM - should probably move this
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
        p.pM = (1-p.m)*eye(h,h) + p.m*repmat(p.K/sum(p.K),h,1);
    end
end

