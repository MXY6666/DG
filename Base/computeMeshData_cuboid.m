% Further classify the leaf faces according to face number. 
% Regard the periodic boudary faces as internal faces.
% Copyright (C) Xiangyi Meng.

function md = computeMeshData_cuboid(msh)

% mesh size
md.meshSize = min(msh.elemLength, [], 'all');

% Classify the internal leaf face IDs
% Based on left and right face number, we can classify the internal leaf
% face IDs into 3 or 27 (without or with local grid refinement) kinds.
isPeriodic = msh.bndTypes == 1;
if any(isPeriodic)
    IndFaceIDs = [msh.intLFaces, msh.bndLFaces{isPeriodic}];
else
    IndFaceIDs = msh.intLFaces;
end
md.nIntLFaces  = length(IndFaceIDs);
md.nIntLFacesx = 0;
md.nIntLFacesy = 0;
md.nIntLFacesz = 0;
md.intLFaces   = cell(27, 3);
temp = cell(1, 30);
for lfn = [2, 4, 6, 11 : 14, 19 : 22, 27 : 30]
    temp{lfn} = IndFaceIDs(msh.faceNums(1, IndFaceIDs) == lfn);
end
for i = 1 : 27
    [lfn, rfn] = getFaceNumber_cuboid(i);
    md.intLFaces{i, 1} = lfn;
    md.intLFaces{i, 2} = rfn;
    if ~isempty(temp{lfn})
        md.intLFaces{i, 3} = temp{lfn}(msh.faceNums(2, temp{lfn}) == rfn);
    end
    nf = length(md.intLFaces{i, 3});
    if any(i == [1, 4 : 11])
        md.nIntLFacesx = md.nIntLFacesx + nf;
    elseif any(i == [2, 12 : 19])
        md.nIntLFacesy = md.nIntLFacesy + nf;
    else
        md.nIntLFacesz = md.nIntLFacesz + nf;
    end
end

% Classify the boundary leaf face IDs
% We further classify each kind of boundary face into 6 kinds according to
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

