function [set1, set2] = getActiveSets(msh, u, Lamb, dU, dLamb, cdt, bs, KKT, newton)

inc = KKT.positivity + KKT.maximum;
u   = reshape(u, [bs.nep, msh.nLElems]);
dU  = reshape(dU, [bs.nb, msh.nLElems]);
% Compute inequality constraints at each limiter point
g = zeros(bs.nep, inc * msh.nLElems);
if (KKT.positivity)
    g(:, 1 : inc : end) = KKT.Umin - u;
end
if (KKT.maximum)
    g(:, KKT.positivity + 1 : inc : end) = u - KKT.Umax;
end  
g = g(:);

% Compute directional derivative of inequality constraints
dgU = zeros(bs.nep, inc * msh.nLElems);
if (KKT.positivity)
    dgU(:, 1 : inc : end) = -cdt * bs.phi{1} * dU; 
end
if (KKT.maximum)
    dgU(:, KKT.positivity + 1 : inc : end) = cdt * bs.phi{1} * dU; 
end
dgU = dgU(:);

% Active sets
alpha_delta = (Lamb >  -g + newton.delta);
beta_delta  = (Lamb <= -g + newton.delta) & (Lamb >= -g - newton.delta);
gamma_delta = (Lamb <  -g - newton.delta);    

func = min(-g, Lamb);
I11beta_delta = beta_delta .* (func >  0) .* (-dgU >  dLamb);
I12beta_delta = beta_delta .* (func >  0) .* (-dgU <= dLamb);
I21beta_delta = beta_delta .* (func <= 0) .* (-dgU >  dLamb);
I22beta_delta = beta_delta .* (func <= 0) .* (-dgU <= dLamb);

set1 = (alpha_delta + I11beta_delta + I22beta_delta) > 0;
set2 = (gamma_delta + I12beta_delta + I21beta_delta) > 0;

end
