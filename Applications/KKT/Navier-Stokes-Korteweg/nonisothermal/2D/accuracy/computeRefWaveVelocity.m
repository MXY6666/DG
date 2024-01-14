function waveVel = computeRefWaveVelocity(msh, md, U, W1, W2, tc, bs)

U  = reshape(U, [bs.nb, 4 * msh.nLElems]);
W1 = reshape(W1, [bs.nb, msh.nLElems]);
W2 = reshape(W2, [bs.nb, msh.nLElems]);
U  = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

waveVel = 0;
for i = 1 : 2 
    faceIDs = md.intLFaces{i, 3};
    leLIDs  = msh.faceElems(1, faceIDs);
    reLIDs  = msh.faceElems(2, faceIDs);
  
    if (i == 1)
        ul  = bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs);
        ur  = bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs); 
        w1l = bs.phi_face{1, 2} * W1(:, leLIDs);
        w1r = bs.phi_face{1, 1} * W1(:, reLIDs);
        w2l = bs.phi_face{1, 2} * W2(:, leLIDs);
        w2r = bs.phi_face{1, 1} * W2(:, reLIDs);        
    elseif (i == 2)
        ul  = bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs);
        ur  = bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * reLIDs); 
        w1l = bs.phi_face{1, 4} * W1(:, leLIDs);
        w1r = bs.phi_face{1, 3} * W1(:, reLIDs);
        w2l = bs.phi_face{1, 4} * W2(:, leLIDs);
        w2r = bs.phi_face{1, 3} * W2(:, reLIDs);
    end
    vl = [ul(:, 2 : 4 : end) ./ ul(:, 1 : 4 : end), ul(:, 3 : 4 : end) ./ ul(:, 1 : 4 : end)];
    vr = [ur(:, 2 : 4 : end) ./ ur(:, 1 : 4 : end), ur(:, 3 : 4 : end) ./ ur(:, 1 : 4 : end)];
    Tl = getTemperature(ul(:, 1 : 4 : end), vl(:, 1 : end / 2), vl(:, end / 2 + 1 : end), ul(:, 4 : 4 : end), w1l, w2l, tc);
    Tr = getTemperature(ur(:, 1 : 4 : end), vr(:, 1 : end / 2), vr(:, end / 2 + 1 : end), ur(:, 4 : 4 : end), w1r, w2r, tc);
    pl = getPressure(ul(:, 1 : 4 : end), Tl);
    pr = getPressure(ur(:, 1 : 4 : end), Tr);

    if (i == 1)
        alphal  = getEigenMax(ul(:, 1 : 4 : end), vl(:, 1 : end / 2), Tl, pl, tc);
        alphar  = getEigenMax(ur(:, 1 : 4 : end), vr(:, 1 : end / 2), Tr, pr, tc);
    elseif (i == 2)
        alphal  = getEigenMax(ul(:, 1 : 4 : end), vl(:, end / 2 + 1 : end), Tl, pl, tc);
        alphar  = getEigenMax(ur(:, 1 : 4 : end), vr(:, end / 2 + 1 : end), Tr, pr, tc);
    end
    waveVel = max(waveVel, max(max(alphal, alphar), [], 'all'));    
end

end