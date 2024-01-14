function [set1, set2] = getActiveSets(rho, p, Lamb, dU, dLamb, r31, c31, v31, KKT, newton)

if (KKT.positivity)
    % Define positivity constraints at each limiter point
    g = KKT.Umin - [rho(:); p(:)];

    % Compute directional derivative of positivity constraints
    M31 = sparse(r31 - 4 * KKT.NDOF, c31, v31, 2 * KKT.LDOF, 4 * KKT.NDOF);
    ndgU = M31 * dU;
    
    % Active sets
    alpha_delta = (Lamb(:) >  -g + newton.delta);
    beta_delta  = (Lamb(:) <= -g + newton.delta) & (Lamb(:) >= -g - newton.delta);
    gamma_delta = (Lamb(:) <  -g - newton.delta);    
    
    func = min(-g, Lamb(:));
    I11beta_delta = beta_delta .* (func >  0) .* (ndgU >  dLamb);
    I12beta_delta = beta_delta .* (func >  0) .* (ndgU <= dLamb);
    I21beta_delta = beta_delta .* (func <= 0) .* (ndgU >  dLamb);
    I22beta_delta = beta_delta .* (func <= 0) .* (ndgU <= dLamb);
    
    set1 = (alpha_delta + I11beta_delta + I22beta_delta) > 0;
    set2 = (gamma_delta + I12beta_delta + I21beta_delta) > 0;
end

end