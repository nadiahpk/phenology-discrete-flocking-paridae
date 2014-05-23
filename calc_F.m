function F = calc_F(p,x,n,y);

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

    % This code is a bit messy still. The basic idea is that I use an adatpive simpson's integration to figure out where good places to put the integration grid points is, and then I run through those points for each arrival distribution selecting the grid with the closest spacing at each time interval.
    x_hi = logninv(1-p.tol,p.mu_g,p.sigma_g);

    % Get a good spacing using the adaptive simpson's
    [intval,b]=adapt_simp(@(x) lognpdf(x,p.mu_g,p.sigma_g),0,x_hi,1e-5); % adequate for most except when we get to mg < 7
    %[intval,b]=adapt_simp(@(x) lognpdf(x,p.mu_g,p.sigma_g),0,x_hi,1e-6); 
    Ncrs = length(b);
    % *** set that 1e-6 to some parameter
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
    F(i,:) = cumtrapz(xk(start_i:end),gf)(end,:);
end
