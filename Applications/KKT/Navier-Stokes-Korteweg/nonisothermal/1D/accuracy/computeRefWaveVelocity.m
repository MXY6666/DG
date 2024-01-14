function waveVel = computeRefWaveVelocity(msh, md, U, W, tc, bs)

U = reshape(U, [bs.nb, 3 * msh.nLElems]);
W = reshape(W, [bs.nb, msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));
% Evaluate density, momentum, energy, derivative of density, temperature 
% and pressure at each face
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
wl = bs.phi_face{1, 2} * W(:, leLIDs);
wr = bs.phi_face{1, 1} * W(:, reLIDs);
vl = ul(2, :) ./ ul(1, :);
vr = ur(2, :) ./ ur(1, :);
Tl = getTemperature(ul(1, :), vl, ul(3, :), wl, tc);
Tr = getTemperature(ur(1, :), vr, ur(3, :), wr, tc);
pl = getPressure(ul(1, :), Tl);
pr = getPressure(ur(1, :), Tr);

alphal  = getEigenMax(ul(1, :), vl, Tl, pl, tc);
alphar  = getEigenMax(ur(1, :), vr, Tr, pr, tc);
waveVel = max(max(alphal, alphar));

end