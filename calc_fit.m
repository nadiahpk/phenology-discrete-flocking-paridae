function fitness = calc_fit(p,x,n,y);

h = length(x);

if nargin == 3
    [W,U] = calc_w(p,x,n);
else
    %W = calc_wdash(p,x,n,y);
    [W,U] = calc_w(p,x,n,y);
end

fitness = max(real(eig(U)));
