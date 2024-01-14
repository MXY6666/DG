% Further classify the leaf faces according to face number. 
% Regard the periodic boudary faces as internal faces.
% Copyright (C) Xiangyi Meng.

function md = computeMeshData_rect(msh)

% mesh size
md.meshSize = min(msh.elemLength, [], 'all');

% Classify the internal leaf faces
% Based on left and right face number, we can classify the internal leaf
% faces into 2 or 10 (without or with local grid refinement) kinds. 
isPeriodic = msh.bndTypes == 1;
if any(isPeriodic)
    IndFaceIDs = [msh.intLFaces, msh.bndLFaces{isPeriodic}];
else
    IndFaceIDs = msh.intLFaces;
end
md.nIntLFaces  = length(IndFaceIDs);
md.nIntLFacesx = 0;
md.nIntLFacesy = 0;
md.intLFaces   = cell(10, 3);
temp = cell(1, 12);
for lfn = [2, 4, 7, 8, 11, 12]
    temp{lfn} = IndFaceIDs(msh.faceNums(1, IndFaceIDs) == lfn);
end
for i = 1 : 10
    [lfn, rfn] = getFaceNumber_rect(i);
    md.intLFaces{i, 1} = lfn;
    md.intLFaces{i, 2} = rfn;
    if ~isempty(temp{lfn})
        md.intLFaces{i, 3} = temp{lfn}(msh.faceNums(2, temp{lfn}) == rfn);
    end
    nf = length(md.intLFaces{i, 3});
    if any(i == [1, 3 : 6])
        md.nIntLFacesx = md.nIntLFacesx + nf;
    else
        md.nIntLFacesy = md.nIntLFacesy + nf;
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


