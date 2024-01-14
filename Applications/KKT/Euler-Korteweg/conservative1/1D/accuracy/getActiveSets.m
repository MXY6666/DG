function [set1, set2] = getActiveSets(msh, rho, Lamb, dU, dLamb, cdt, bs, KKT, newton)

if (KKT.positivity)
    % Define positivity constraints at each limiter point
    g = KKT.Umin - rho(:);

    % Compute directional derivative of positivity constraints
    dgU = -cdt * bs.phi{1} * reshape(dU(1 : KKT.NDOF), [bs.nb, msh.nLElems]);
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

end
