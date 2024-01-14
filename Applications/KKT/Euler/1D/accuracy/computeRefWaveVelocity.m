function waveVel = computeRefWaveVelocity(msh, md, U, tc, bs)

U = reshape(U, [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));
% Evaluate density, momentum, energy and pressure at each face
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = getPressure(ur(1, :), ur(2, :), ur(3, :), tc);

alphal  = getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc);
alphar  = getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc);
waveVel = max(max(alphal, alphar));

end