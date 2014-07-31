function F = calc_F(p,x,n,y);

% -- F = calc_F(p,x,n,y);
% -- F = calc_F(p,x,n);
%
%
% The purpose of this function is to calculate the
% probability that a juvenile from each habitat will be
% first in the territory queue in each habitat.
%
%
% INPUTS
%
% p: The dictionary of parameter values. See marsh2.m for an
% example of how to specify these.
%
% x: The hatching-date strategy being pursued by the
% population ("resident strategy"). Row vector.
%
% n: The population size. Row vector.
%
% y: The hatching-date strategy being pursued by the
% variant ("mutant strategy"). Row vector. Optional.
%
%
% OUTPUTS
%
% F: The probability of being first in the queue for an
% individual travelling from habitat i (row) to habitat j
% (column). Depending upon whether y is specified, it is
% calculated either for a resident or mutant invader.
%


h = length(x); % Number of habitats, traits

intsteps_method = 1; % 0: simple evenly spaced
                     % 1: using the adaptive simpson's with varying size

if intsteps_method == 0 % Evenly spaced option

    Ncrs = p.Ncrs;
    tol = p.tol;
    if nargin == 3
        x_hi = max(x+logninv(1-p.tol,p.mu_g,p.sigma_g));
        x_lo = min(x);
    else
        x_hi = max([x,y]+logninv(1-p.tol,p.mu_g,p.sigma_g));
        x_lo = min([x,y]);
    end

    xk = linspace(x_lo,x_hi,Ncrs);

elseif intsteps_method == 1 % More complicated option

    % This more-complicated approach became necessary when
    % the skew of the arrival distribution was very strong
    % and the peak very sharp.
    %
    % Here an adaptive simpson's integration is used to
    % figure out where are good places to put the
    % integration grid points. The code then runs through
    % those points for each arrival distribution, selecting
    % the grid with the closest spacing at each time interval.

    x_hi = logninv(1-p.tol,p.mu_g,p.sigma_g);

    % Get a good spacing using the adaptive simpson's
    [intval,b]=adapt_simp(@(x) lognpdf(x,p.mu_g,p.sigma_g),0,x_hi,p.adaptsimp);
    Ncrs = length(b);

    % Determine spacing
    bb = b; bb(end)=[]; bb = [NaN,bb]; b_spacing = bb-b; 

    % Generate a grid for each of our arrival distributions
    if nargin == 4
        xx = [x,y];
        xx = repmat([x,y]',1,Ncrs)+repmat(b,h*2,1);
        nrows = h*2;
    else
        xx = repmat(x',1,Ncrs)+repmat(b,h,1);
        nrows = h;
    end

    % Now go through selecting the most closely-spaced grid
    % across the entire interval of interest
    col = Ncrs;
    [xk(1),row] = max(xx(:,end));
    xx(row,col) = NaN;
    while not(isempty(find(not(isnan(xx)))));

        [vals,col2] = max(xx');
        nextspacing = b_spacing(col2);
        spacing = nextspacing(row);

        row2 = 1:nrows;
        remind = [find(vals <= vals(row)),find(isnan(vals))];
        col2(remind) = [];
        row2(remind) = [];
        if isempty(row2); % If next in row is next in order
            col = col-1;
        else
            nextspacing = b_spacing(col2);
            % Retain only those rows that are more closely spaced than the current spacing
            remind = [find(nextspacing<=spacing),find(isnan(nextspacing))];
            row2(remind) = [];
            if isempty(row2)
                col = col-1;
            else
                % Choose the one with the closest spacing
                col2(remind) = [];
                nextspacing(remind) = [];
                [a,ind] = max(nextspacing);
                row = row2(ind);
                col = col2(ind);
            end
        end
        xk = [xx(row,col),xk];
        xx(find(xx >= xx(row,col))) = NaN;
    end

    Ncrs = length(xk); % Update length of course grid
end

% Find f,g on the grid
[f,gg]=calc_fg(p,x,n,xk);

if nargin == 4
    % Recalculate gg for our mutant because we're finding F'
    gg = lognpdf(repmat(xk',1,h)-repmat(y,Ncrs,1),p.mu_g,p.sigma_g);
end

% Integrate from min(x) to the point at which gg is close enough to zero
F = zeros(h,h);
for i = 1:h
    if nargin == 4
        [a,start_i] = min(abs(xk-y(i))); 
    else
        [a,start_i] = min(abs(xk-x(i))); 
    end
    gf = repmat(gg(start_i:end,i),1,h).*f(start_i:end,:);
    F(i,:) = cumtrapz(xk(start_i:end)',gf)(end,:);
end
