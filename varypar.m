function varypar(p,par_name,parV,par_col,mm0,n0);

% -- varypar(p,par_name,parV,par_col,mm0,n0)
%
% The purpose of this function is to automate the
% exploration of the effects of each of the parameters upon
% the ESS. It will produce a .dat file with the name:
% 'output'par_name'_'par_col'.dat'
%
% Note that I have hardcoded in the variation of the optimal
% hatching time to match the figures in the paper.
%
% Example of use, explore effect of K in early habitat:
% Get parameter-value dictionary p
%   default_parameter_values 
% Values of K at which we'd like to evaluate the ESS
%   parV = [linspace(20,80,13)',linspace(80,20,13)'] 
% Initial guess of the mismatch 
%   mm0=[-8.1471e+00     -8.0800e+00] 
% Initial guess of the population sizes
%   n0=[1.8549e+01      7.4211e+01];
%   varypar(p,'K',parV,1,mm0,n0)
%
%
% INPUTS 
%
% p: The dictionary of parameter values
%
% par_name: The parameter name as a string. For example, if
% we want to vary p.K, we'd enter have here par_name = 'K'
%
% parV: A vector or matrix of the parameter values at which
% we will evaluate the ESS. Each row is taken to be another
% parameter value, so e.g. if we have a parameter like K
% that is two values (i.e. [K in early habitat, K in late
% habitat]), then this would be an n x 2 matrix.
%
% par_col: Indicates which column of the parV is varying. My
% convention has been to set the earlier columns to the
% earlier habitats, so e.g. if we are interested in sigma in the 
% early habitat (as opposed to the late habitat) we'd say par_col = 1 
%
% mm0: An initial guess for the mismatch in each habitat at
% our start value
%
% n0: An initial guess for the population size in each habitat at
% our start value
%

par_n = size(parV,1);

% I've chosen this to reproduce Figure 1, however you may
% wish to change this range
x1V = [linspace(171,161,11),linspace(159,131,15)];
x1_n = 26;

% Get start value
eval(['p.',par_name,' = [',num2str(parV(1,:)),'];']);
eval(['p.x_opt(1) = ',num2str(x1V(1)),';']);
p = updated_pm(p); % Updates after change to p.K

stormm1 = zeros(x1_n,par_n);
stormm2 = stormm1;

str = ['output',par_name,'_',num2str(par_col),'.dat'];
fid = fopen(str,'w');

fprintf(fid,'# Generated with varypar.m\n');
fprintf(fid,'#\n');
fprintf(fid,'# Default parameter values were:\n');
for [val,key] = p;
    if size(val,1) == 1; % Don't print pM basically
        str = ['# p.',key,' = ',num2str(val),';\n'];
        fprintf(fid,str);
    end
end
fprintf(fid,'#\n');
fprintf(fid,'# Column entries are:\n');
fprintf(fid,'# 1:x_1 \t 2:parval \t 3:x_1-x_1^* \t 4:x_2-x_2^* \t 5:n1 \t 6:n2 \t 7:flow11 \t 8:flow21 \t 9:flow12 \t 10:flow22 \n');
fprintf(fid,'\n');

for ind_x1 = 1:x1_n;

    x1 = x1V(ind_x1);
    p.x_opt(1) = x1;

    x = p.x_opt+mm0;
    n = n0;

    for ind_par = 1:par_n

        par = parV(ind_par,:); 
        eval(['p.',par_name,' = par;']);

        p = updated_pm(p); % Updates after change to p.K
        x = calc_x(p,x,n);
        n = calc_n(p,x);
        %x = calc_x(p,x,p.K);

        % Store mismatches
        mm = x-p.x_opt;
        stormm1(ind_x1,ind_par) = mm(1);
        stormm2(ind_x1,ind_par) = mm(2);
        
        flow = calc_flow(p,x,n);
        flow = reshape(flow,1,4);

        fprintf(fid,'%0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \n',[x1,par(par_col),mm,n,flow]);
        fflush(fid);

        % Stor for next x1 step
        disp(['x1 = ',num2str(x1),', par = ',num2str(par),', x = ',num2str(x),'. ']);
        fflush(stdout);
        if ind_par == 1
            mm0 = mm;
            n0 = n;
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
