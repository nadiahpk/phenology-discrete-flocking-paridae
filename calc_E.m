function E = calc_E(p,xC)

% x is a vector, indices correspond to habitat
% x_opt is a vector, indices correspond to habitat

E = p.E_0.*exp(-(xC-p.x_opt).^2./(2*p.sigma_E.^2));

