clear all;
default_parameter_values;

% Read singular strategies from file into M
fid = fopen('outputK_1.dat');
fskipl(fid, 22); % Comments lines
str = '%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g';
[M, nM] = fscanf(fid,str,[10, Inf]);
M=M';
M = M(:,1:6);
fclose(fid);
rows = size(M)(1);

parV = [M(:,2),100-M(:,2)];

% Open a storage file and write comments
outfile_name = ['checkstabK_1.dat'];
fid = fopen(outfile_name,'w');
fprintf(fid,'# Generated with scriptcheckstab.m\n');
fprintf(fid,['# for data from file: outputK_1.dat\n']);
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
fprintf(fid,'# 1:x_1 \t 2:parval \t 3:x_1-x_1^* \t 4:x_2-x_2^* \t 5:n1 \t 6:n2 \t 7:eigH \t 8:eigJ \t 9:eigJs \n');
fprintf(fid,'\n');

start_ind = 1; % -- MODIFY THIS if it gets stuck somewhere
x1_prev = M(start_ind,1); % For checking when we need to insert a break into the dat file

% Go through each row checking stability
for row = start_ind:rows
    
    % Get steady state trait and population
    p.x_opt(1) = M(row,1);
    x = p.x_opt + M(row,3:4);
    n = M(row,5:6);

    % Get the parameter values for this steady state
    par_val = parV(row,:);
    eval(['p.K = [',num2str(par_val),'];']);
    p = updated_pm(p); % Needed after changes to p.K

    % Check evolutionary stability
    [eigH,eigJ,eigJs,Hess,Jac]=check_stab(p,x,n);
    % debugging eigH = rand; eigJ = rand; eigJs = rand; Hess = rand; Jac = rand;

    % Write evolutionary stability results to file:
    % Insert a break for change in x1
    if p.x_opt(1) ~= x1_prev 
        fprintf(fid,'\n');
        x1_prev = p.x_opt(1);
    end
    % Print line of data
    fprintf(fid,'%0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \t %0.4e \n',[M(row,:),eigH,eigJ,eigJs]);
    fflush(fid);

    % Display to screen so I can watch progress
    disp(['x(1) = ',num2str(p.x_opt(1)),', par = ',num2str(par_val),', eigH = ',num2str(eigH),', eigJ = ',num2str(eigJ),', eigJs = ',num2str(eigJs),'. ']);
    fflush(stdout);

end

fclose(fid);

