function waveVel = computeRefWaveVelocity(msh, md, U, tc, bs)

U = reshape(U, [bs.nb, 4 * msh.nLElems]);
U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

waveVel = 0;
for i = 1 : 2 
    faceIDs = md.intLFaces{i, 3};
    leLIDs  = msh.faceElems(1, faceIDs);
    reLIDs  = msh.faceElems(2, faceIDs);

    if (i == 1)
        ul = bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs);
        ur = bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs); 
        vl = ul(:, 2 : 4 : end) ./ ul(:, 1 : 4 : end);
        vr = ur(:, 2 : 4 : end) ./ ur(:, 1 : 4 : end);
    elseif (i == 2)
        ul = bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs);
        ur = bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * reLIDs); 
        vl = ul(:, 3 : 4 : end) ./ ul(:, 1 : 4 : end);
        vr = ur(:, 3 : 4 : end) ./ ur(:, 1 : 4 : end);
    end
    pl = getPressure(ul(:, 1 : 4 : end), ul(:, 2 : 4 : end), ul(:, 3 : 4 : end), ul(:, 4 : 4 : end), tc);
    pr = getPressure(ur(:, 1 : 4 : end), ur(:, 2 : 4 : end), ur(:, 3 : 4 : end), ur(:, 4 : 4 : end), tc);
    
    alphal  = getEigenMax(ul(:, 1 : 4 : end), vl, pl, tc);
    alphar  = getEigenMax(ur(:, 1 : 4 : end), vr, pr, tc);
    waveVel = max(waveVel, max(max(alphal, alphar), [], 'all'));    
end

end