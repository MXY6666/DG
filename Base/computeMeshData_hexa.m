% Further classify the leaf faces according to face number. 
% Regard the periodic boudary faces as internal faces.
% Copyright (C) Xiangyi Meng.

function md = computeMeshData_hexa(msh)

% mesh size
md.meshSize = min(msh.elemDiam);

% Classify the internal leaf face IDs
% Based on left and right face number, we can classify the internal leaf
% face IDs into 36 or 324 (without or with local grid refinement) kinds.
isPeriodic = msh.bndTypes == 1;
if any(isPeriodic)
    IndFaceIDs = [msh.intLFaces, msh.bndLFaces{isPeriodic}];
else
    IndFaceIDs = msh.intLFaces;
end
md.nIntLFaces = length(IndFaceIDs);
md.intLFaces  = cell(324, 3, 8);
md.nifs = zeros(1, 8);
temp = cell(1, 30);
for f2f = 1 : 8
    faceIDs = IndFaceIDs(msh.faceR2L(IndFaceIDs) == f2f);
    if ~isempty(faceIDs)
        md.nifs(f2f) = length(faceIDs);
        for lfn = 1 : 30
            temp{lfn} = faceIDs(msh.faceNums(1, faceIDs) == lfn);
        end
        for i = 1 : 324
            if (i <= 36)
                lfn = ceil(i / 6); rfn = i - 6 * (lfn - 1);
            else
                if (i <= 180)
                    j = i - 36;
                    lfn = ceil(j / 24); rfn = j - 24 * (lfn - 1) + 6;
                else
                    j = i - 180;
                    lfn = 6 + ceil(j / 6); rfn = j - 6 * (lfn - 7);
                end
            end  
            md.intLFaces{i, 1, f2f} = lfn;
            md.intLFaces{i, 2, f2f} = rfn;
            if ~isempty(temp{lfn})
                md.intLFaces{i, 3, f2f} = temp{lfn}(msh.faceNums(2, temp{lfn}) == rfn);
            end
        end        
    end
end

% Classify the boundary leaf face IDs
% We further classify each kind of boundary face into 4 kinds according to
% left face number. 
nBndTypes     = length(msh.bndTypes);
md.nBndLFaces = msh.nBndLFaces;
md.bndLFaces  = cell(6, nBndTypes);
for i = 1 : nBndTypes
    if (msh.bndTypes(i) == 1)
        md.nBndLFaces(i) = 0;
    else
        for lfn = 1 : 6
            md.bndLFaces{lfn, i} = msh.bndLFaces{i}(msh.faceNums(1, msh.bndLFaces{i}) == lfn);
        end
    end
end

end








