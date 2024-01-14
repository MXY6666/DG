% Further classify the leaf faces according to face number. 
% Regard the periodic boudary faces as internal faces. 
% Copyright (C) Xiangyi Meng.

function md = computeMeshData_quad(msh)

% mesh size
md.meshSize = min(msh.elemDiam);

% Classify the internal leaf face IDs
% Based on left and right face number, we can classify the internal leaf
% face IDs into 16 or 80 (without or with local grid refinement) kinds. 
isPeriodic = msh.bndTypes == 1;
if any(isPeriodic)
    IndFaceIDs = [msh.intLFaces, msh.bndLFaces{isPeriodic}];
else
    IndFaceIDs = msh.intLFaces;
end
md.nIntLFaces = length(IndFaceIDs);
md.intLFaces  = cell(80, 3);
temp = cell(1, 12);
for lfn = 1 : 12
    temp{lfn} = IndFaceIDs(msh.faceNums(1, IndFaceIDs) == lfn);
end
for i = 1 : 80
    if (i <= 16)
        lfn = ceil(i / 4); rfn = i - 4 * (lfn - 1);
    else
        if (i <= 48)
            j = i - 16;
            lfn = ceil(j / 8); rfn = j - 8 * (lfn - 1) + 4;
        else
            j = i - 48;
            lfn = 4 + ceil(j / 4); rfn = j - 4 * (lfn - 5);
        end
    end  
    md.intLFaces{i, 1} = lfn;
    md.intLFaces{i, 2} = rfn;
    if ~isempty(temp{lfn})
        md.intLFaces{i, 3} = temp{lfn}(msh.faceNums(2, temp{lfn}) == rfn);
    end
end

% Classify the boundary leaf face IDs
% We further classify each kind of boundary face into 4 kinds according to
% left face number. 
nBndTypes     = length(msh.bndTypes);
md.nBndLFaces = msh.nBndLFaces;
md.bndLFaces  = cell(4, nBndTypes);
for i = 1 : nBndTypes
    if (msh.bndTypes(i) == 1)
        md.nBndLFaces(i) = 0;
    else
        for lfn = 1 : 4
            md.bndLFaces{lfn, i} = msh.bndLFaces{i}(msh.faceNums(1, msh.bndLFaces{i}) == lfn);
        end
    end
end

end



