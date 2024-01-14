% Suppose the level difference between neighboring elements is no more than
% 1 when refining and coarsening
% Copyright (C) Xiangyi Meng.

function [msh, U] = refineAndCoarsenCuboidMesh(msh, candReIDs, candCeIDs, nv, U, c2f, f2c, layout)

if (nargin < 1)
    error('Not enough arguments')
end

if (msh.type ~= 302) 
    error('The structure of the given mesh is not adaptive')
end

if (nargin < 2) 
    candReIDs = [];
end

if (nargin < 3) 
    candCeIDs = [];
end
if isempty(candReIDs) && isempty(candCeIDs)
    return
end

if (nargin < 4) || isempty(nv)
    nv = 1;
end

if (nargin < 5) || isempty(U)
    U = zeros(4, msh.nLElems);
end
if (size(U, 2) ~= msh.nLElems * nv)
    error('The given mesh and the solution coefficients does not coinside')
end

if (nargin < 6) || isempty(c2f)
    c2f = getC2FMatrix_cuboid(1);
end

if (nargin < 7) || isempty(f2c)
    f2c = getF2CMatrix_cuboid(1);
end

if (nargin < 8) || isempty(layout)
    layout = 2;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

%**************************************************************************
% Firstly refine the elements. Note that we may add more elements to refine
% to ensure the level difference between neighboring elements is no more
% than 1
%**************************************************************************
% Delete the candidate elements whose level is greater or equal to maxLevel
reIDs = candReIDs(msh.elemLevel(candReIDs) < msh.maxLevel);
nre = length(reIDs);
shouldRefine = false(1, msh.nElems);
shouldRefine(reIDs) = true;

% Make sure the level difference of adjacent elements is no more than 1 
% by increasing the number of refined elements
reIDs(nre + 1 : floor(msh.nLElems / 3)) = 0; % preallocation
i = 1;
while i <= nre
    % fields of the current loop element
    ID = reIDs(i); 
    level     = msh.elemLevel(ID);
    elemFaces = msh.elemFaces(:, ID);   
    faceTypes = msh.faceType(elemFaces);
    faceElems = msh.faceElems(:, elemFaces);
    
    % Determine if refine the first adjacent element
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        ID1 = faceElems(1, 1);
        if (msh.elemLevel(ID1) < level) && ~shouldRefine(ID1)
            shouldRefine(ID1) = true;
            nre = nre + 1;
            reIDs(nre) = ID1;
        end       
    end
    % Determine if refine the second adjacent element
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        ID2 = faceElems(2, 2);
        if (msh.elemLevel(ID2) < level) && ~shouldRefine(ID2)
            shouldRefine(ID2) = true;
            nre = nre + 1;
            reIDs(nre) = ID2;
        end         
    end
    % Determine if refine the third adjacent element
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        ID3 = faceElems(1, 3);
        if (msh.elemLevel(ID3) < level) && ~shouldRefine(ID3)
            shouldRefine(ID3) = true;
            nre = nre + 1;
            reIDs(nre) = ID3;
        end         
    end  
    % Determine if refine the fourth adjacent element
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        ID4 = faceElems(2, 4);
        if (msh.elemLevel(ID4) < level) && ~shouldRefine(ID4)
            shouldRefine(ID4) = true;
            nre = nre + 1;
            reIDs(nre) = ID4;
        end         
    end    
    % Determine if refine the fifth adjacent element
    if (faceTypes(5) == 0 || faceTypes(5) == 1)
        ID5 = faceElems(1, 5);
        if (msh.elemLevel(ID5) < level) && ~shouldRefine(ID5)
            shouldRefine(ID5) = true;
            nre = nre + 1;
            reIDs(nre) = ID5;
        end         
    end    
    % Determine if refine the sixth adjacent element
    if (faceTypes(6) == 0 || faceTypes(6) == 1)
        ID6 = faceElems(2, 6);
        if (msh.elemLevel(ID6) < level) && ~shouldRefine(ID6)
            shouldRefine(ID6) = true;
            nre = nre + 1;
            reIDs(nre) = ID6;
        end         
    end
    i = i + 1;
end
reIDs(nre + 1 : end) = [];

% We should first refine the elements of level 0, then level 1, then ...
% finally level (maxLvel - 1)
reLevels = msh.elemLevel(reIDs);
sreIDs = zeros(1, nre);
s = 0;
for i = 0 : msh.maxLevel - 1
    temp = reIDs(reLevels == i);
    l = length(temp);
    sreIDs(s + 1 : s + l) = temp;
    s = s + l;
end
if (s ~= nre)
    error('Wrong size after sorting')
end
reIDs = sreIDs;

% Then we refine the elements by adding new elements and faces
% preallocation for elements
msh.elemCenter(:, msh.nElems + 1 : msh.nElems + 8 * nre)   = 0;
msh.elemLength(:, msh.nElems + 1 : msh.nElems + 8 * nre)   = 0;
msh.elemSize(msh.nElems + 1 : msh.nElems + 8 * nre)        = 0;
msh.elemFaces(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0;
msh.elemLevel(msh.nElems + 1 : msh.nElems + 8 * nre)       = 0;
msh.elemParent(msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemChildren(:, msh.nElems + 1 : msh.nElems + 8 * nre) = 0;
msh.elemJac(:, msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemJxix(:, msh.nElems + 1 : msh.nElems + 8 * nre)     = 0;
msh.elemJetay(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0;
msh.elemJzetaz(:, msh.nElems + 1 : msh.nElems + 8 * nre)   = 0;
% preallocation for faces
msh.faceNormalx(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)  = 0; 
msh.faceNormaly(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)  = 0; 
msh.faceNormalz(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)  = 0; 
msh.faceSize(msh.nFaces + 1 : msh.nFaces + 36 * nre)        = 0; 
msh.faceType(msh.nFaces + 1 : msh.nFaces + 36 * nre)        = 0; 
msh.faceElems(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)    = 0; 
msh.faceNums(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)     = 0; 
msh.faceLevel(msh.nFaces + 1 : msh.nFaces + 36 * nre)       = 0; 
msh.faceParent(msh.nFaces + 1 : msh.nFaces + 36 * nre)      = 0; 
msh.faceChildren(:, msh.nFaces + 1 : msh.nFaces + 36 * nre) = 0; 
msh.faceJac(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)      = 0; 
% preallocation for solution
UT = zeros(size(U, 1), (msh.nElems + 8 * nre) * nv);
if (layout == 1)
    U = U(:, (0 : nv - 1)' * msh.nLElems + (1 : msh.nLElems));
end
UT(:, (-nv + 1 : 0)' + nv * msh.LElems) = U;

faceIDs1 = zeros(3, 2, 2);
faceIDs2 = zeros(2, 3, 2);
faceIDs3 = zeros(2, 2, 3);
for i = 1 : nre
    % current loop element and its faces
    ID = reIDs(i);
    ct = msh.elemCenter(:, ID);
    h  = msh.elemLength(:, ID);
    si = msh.elemSize(ID);
    elemFaces    = msh.elemFaces(:, ID);
    faceNormalxs = msh.faceNormalx(elemFaces);
    faceNormalys = msh.faceNormaly(elemFaces);
    faceNormalzs = msh.faceNormalz(elemFaces);
    faceTypes    = msh.faceType(elemFaces);
    faceElems    = msh.faceElems(:, elemFaces);
    faceLevels   = msh.faceLevel(elemFaces);
    faceChildren = msh.faceChildren(:, elemFaces);    
        
    % creation of subelements
    elemChildren = reshape(msh.nElems + 1 : msh.nElems + 8, [2, 2, 2]);
    msh.elemChildren(:, ID)         = elemChildren(:); 
    msh.elemCenter(:, elemChildren) = [ct(1) - h(1) / 4, ct(1) + h(1) / 4, ct(1) - h(1) / 4, ct(1) + h(1) / 4, ct(1) - h(1) / 4, ct(1) + h(1) / 4, ct(1) - h(1) / 4, ct(1) + h(1) / 4; ...
                                       ct(2) - h(2) / 4, ct(2) - h(2) / 4, ct(2) + h(2) / 4, ct(2) + h(2) / 4, ct(2) - h(2) / 4, ct(2) - h(2) / 4, ct(2) + h(2) / 4, ct(2) + h(2) / 4; ...
                                       ct(3) - h(3) / 4, ct(3) - h(3) / 4, ct(3) - h(3) / 4, ct(3) - h(3) / 4, ct(3) + h(3) / 4, ct(3) + h(3) / 4, ct(3) + h(3) / 4, ct(3) + h(3) / 4]; 
    msh.elemLength(:, elemChildren) = repmat(h / 2, [1, 8]);
    msh.elemSize(elemChildren)      = si / 8; 
    msh.elemLevel(elemChildren)     = msh.elemLevel(ID) + 1;
    msh.elemParent(elemChildren)    = ID;
    msh.elemJac(elemChildren)       = si / 64;
    msh.elemJxix(elemChildren)      = h(2) * h(3) / 16;
    msh.elemJetay(elemChildren)     = h(1) * h(3) / 16;
    msh.elemJzetaz(elemChildren)    = h(1) * h(2) / 16;
        
    % creation of subfaces
    % first face
    if (faceChildren(1, 1) == 0)
        faceChildren(:, 1) = (msh.nFaces + 1 : msh.nFaces + 4)';
        msh.faceChildren(:, elemFaces(1)) = faceChildren(:, 1); 
        
        faceIDs1(1, :, :) = reshape(faceChildren(:, 1), [1, 2, 2]);
        if (faceTypes(1) == 0 || faceTypes(1) == 1)
            msh.faceElems(1, faceChildren(:, 1)) = faceElems(1, 1);
            msh.faceElems(2, faceChildren(:, 1)) = reshape(elemChildren(1, :, :), [1, 4]);
            msh.faceNums( 1, faceChildren(:, 1)) = 11 : 14;
            msh.faceNums( 2, faceChildren(:, 1)) = 1;            
        else
            msh.faceElems(1, faceChildren(:, 1)) = reshape(elemChildren(1, :, :), [1, 4]);
            msh.faceNums( 1, faceChildren(:, 1)) = 1; 
        end
        msh.faceNormalx(faceChildren(:, 1)) = faceNormalxs(1);
        msh.faceSize(faceChildren(:, 1))    = h(2) * h(3) / 4;
        msh.faceType(faceChildren(:, 1))    = faceTypes(1);
        msh.faceLevel(faceChildren(:, 1))   = faceLevels(1) + 1;
        msh.faceParent(faceChildren(:, 1))  = elemFaces(1);
        msh.faceJac(faceChildren(:, 1))     = h(2) * h(3) / 16;
        msh.nFaces = msh.nFaces + 4;        
    else
        faceIDs1(1, :, :) = reshape(faceChildren(:, 1), [1, 2, 2]);
        % It must be an internal face in this case
        msh.faceElems(2, faceChildren(:, 1)) = reshape(elemChildren(1, :, :), [1, 4]);
        msh.faceNums( 2, faceChildren(:, 1)) = 1; 
    end
    % second face
    if (faceChildren(1, 2) == 0)
        faceChildren(:, 2) = (msh.nFaces + 1 : msh.nFaces + 4)';
        msh.faceChildren(:, elemFaces(2)) = faceChildren(:, 2);    
 
        faceIDs1(3, :, :) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [1, 2, 2]); 
        msh.faceElems(1, faceChildren(:, 2)) = reshape(elemChildren(2, :, :), [1, 4]);
        msh.faceNums( 1, faceChildren(:, 2)) = 2;
        if (faceTypes(2) == 0 || faceTypes(2) == 1) 
            msh.faceElems(2, faceChildren(:, 2)) = faceElems(2, 2);           
            msh.faceNums( 2, faceChildren(:, 2)) = 7 : 10;
        end  
        msh.faceNormalx(faceChildren(:, 2)) = faceNormalxs(2);
        msh.faceSize(faceChildren(:, 2))    = h(2) * h(3) / 4;
        msh.faceType(faceChildren(:, 2))    = faceTypes(2);
        msh.faceLevel(faceChildren(:, 2))   = faceLevels(2) + 1;
        msh.faceParent(faceChildren(:, 2))  = elemFaces(2);
        msh.faceJac(faceChildren(:, 2))     = h(2) * h(3) / 16;
        msh.nFaces = msh.nFaces + 4;        
    else
        faceIDs1(3, :, :) = reshape(faceChildren(:, 2), [1, 2, 2]);
        % It must be an internal face in this case
        msh.faceElems(1, faceChildren(:, 2)) = reshape(elemChildren(2, :, :), [1, 4]);
        msh.faceNums( 1, faceChildren(:, 2)) = 2;         
    end   
    % third face
    if (faceChildren(1, 3) == 0)
        faceChildren(:, 3) = (msh.nFaces + 1 : msh.nFaces + 4)';
        msh.faceChildren(:, elemFaces(3)) = faceChildren(:, 3);  
        
        faceIDs2(:, 1, :) = reshape(msh.nFaces + [1, 3, 2, 4], [2, 1, 2]); 
        if (faceTypes(3) == 0 || faceTypes(3) == 1)
            msh.faceElems(1, faceChildren(:, 3)) = faceElems(1, 3);
            msh.faceElems(2, faceChildren(:, 3)) = elemChildren([1, 5, 2, 6]);
            msh.faceNums( 1, faceChildren(:, 3)) = 19 : 22;
            msh.faceNums( 2, faceChildren(:, 3)) = 3; 
        else
            msh.faceElems(1, faceChildren(:, 3)) = elemChildren([1, 5, 2, 6]);
            msh.faceNums( 1, faceChildren(:, 3)) = 3; 
        end 
        msh.faceNormaly(faceChildren(:, 3)) = faceNormalys(3);
        msh.faceSize(faceChildren(:, 3))    = h(1) * h(3) / 4;
        msh.faceType(faceChildren(:, 3))    = faceTypes(3);
        msh.faceLevel(faceChildren(:, 3))   = faceLevels(3) + 1;
        msh.faceParent(faceChildren(:, 3))  = elemFaces(3);
        msh.faceJac(faceChildren(:, 3))     = h(1) * h(3) / 16;
        msh.nFaces = msh.nFaces + 4;        
    else
        faceIDs2(:, 1, :) = reshape(faceChildren([1, 3, 2, 4], 3), [2, 1, 2]); 
        % It must be an internal face in this case
        msh.faceElems(2, faceChildren(:, 3)) = elemChildren([1, 5, 2, 6]);
        msh.faceNums( 2, faceChildren(:, 3)) = 3;         
    end    
    % fourth face
    if (faceChildren(1, 4) == 0)
        faceChildren(:, 4) = (msh.nFaces + 1 : msh.nFaces + 4)';
        msh.faceChildren(:, elemFaces(4)) = faceChildren(:, 4);    
        
        faceIDs2(:, 3, :) = reshape(msh.nFaces + [1, 3, 2, 4], [2, 1, 2]);
        msh.faceElems(1, faceChildren(:, 4)) = elemChildren([3, 7, 4, 8]);
        msh.faceNums( 1, faceChildren(:, 4)) = 4;        
        if (faceTypes(4) == 0 || faceTypes(4) == 1)
            msh.faceElems(2, faceChildren(:, 4)) = faceElems(2, 4);
            msh.faceNums( 2, faceChildren(:, 4)) = 15 : 18;
        end 
        msh.faceNormaly(faceChildren(:, 4)) = faceNormalys(4);
        msh.faceSize(faceChildren(:, 4))    = h(1) * h(3) / 4;
        msh.faceType(faceChildren(:, 4))    = faceTypes(4);
        msh.faceLevel(faceChildren(:, 4))   = faceLevels(4) + 1;
        msh.faceParent(faceChildren(:, 4))  = elemFaces(4);
        msh.faceJac(faceChildren(:, 4))     = h(1) * h(3) / 16;
        msh.nFaces = msh.nFaces + 4;         
    else
        faceIDs2(:, 3, :) = reshape(faceChildren([1, 3, 2, 4], 4), [2, 1, 2]);
        % It must be an internal face in this case
        msh.faceElems(1, faceChildren(:, 4)) = elemChildren([3, 7, 4, 8]);
        msh.faceNums( 1, faceChildren(:, 4)) = 4;         
    end   
    % fifth face
    if (faceChildren(1, 5) == 0)
        faceChildren(:, 5) = (msh.nFaces + 1 : msh.nFaces + 4)';
        msh.faceChildren(:, elemFaces(5)) = faceChildren(:, 5); 
        
        faceIDs3(:, :, 1) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [2, 2]);
        if (faceTypes(5) == 0 || faceTypes(5) == 1)
            msh.faceElems(1, faceChildren(:, 5)) = faceElems(1, 5);
            msh.faceElems(2, faceChildren(:, 5)) = reshape(elemChildren(:, :, 1), [1, 4]);
            msh.faceNums( 1, faceChildren(:, 5)) = 27 : 30;
            msh.faceNums( 2, faceChildren(:, 5)) = 5; 
        else
            msh.faceElems(1, faceChildren(:, 5)) = reshape(elemChildren(:, :, 1), [1, 4]);
            msh.faceNums( 1, faceChildren(:, 5)) = 5; 
        end 
        msh.faceNormalz(faceChildren(:, 5)) = faceNormalzs(5);
        msh.faceSize(faceChildren(:, 5))    = h(1) * h(2) / 4;
        msh.faceType(faceChildren(:, 5))    = faceTypes(5);
        msh.faceLevel(faceChildren(:, 5))   = faceLevels(5) + 1;
        msh.faceParent(faceChildren(:, 5))  = elemFaces(5);
        msh.faceJac(faceChildren(:, 5))     = h(1) * h(2) / 16;
        msh.nFaces = msh.nFaces + 4;        
    else
        faceIDs3(:, :, 1) = reshape(faceChildren(:, 5), [2, 2]);
        % It must be an internal face in this case
        msh.faceElems(2, faceChildren(:, 5)) = reshape(elemChildren(:, :, 1), [1, 4]);
        msh.faceNums( 2, faceChildren(:, 5)) = 5;          
    end    
    % sixth face
    if (faceChildren(1, 6) == 0)
        faceChildren(:, 6) = (msh.nFaces + 1 : msh.nFaces + 4)';
        msh.faceChildren(:, elemFaces(6)) = faceChildren(:, 6);
        
        faceIDs3(:, :, 3) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [2, 2]); 
        msh.faceElems(1, faceChildren(:, 6)) = reshape(elemChildren(:, :, 2), [1, 4]);
        msh.faceNums( 1, faceChildren(:, 6)) = 6;        
        if (faceTypes(6) == 0 || faceTypes(6) == 1)
            msh.faceElems(2, faceChildren(:, 6)) = faceElems(2, 6);
            msh.faceNums( 2, faceChildren(:, 6)) = 23 : 26;            
        end  
        msh.faceNormalz(faceChildren(:, 6)) = faceNormalzs(6);
        msh.faceSize(faceChildren(:, 6))    = h(1) * h(2) / 4;
        msh.faceType(faceChildren(:, 6))    = faceTypes(6);
        msh.faceLevel(faceChildren(:, 6))   = faceLevels(6) + 1;
        msh.faceParent(faceChildren(:, 6))  = elemFaces(6);
        msh.faceJac(faceChildren(:, 6))     = h(1) * h(2) / 16;
        msh.nFaces = msh.nFaces + 4;         
    else
        faceIDs3(:, :, 3) = reshape(faceChildren(:, 6), [2, 2]);
        % It must be an internal face in this case
        msh.faceElems(1, faceChildren(:, 6)) = reshape(elemChildren(:, :, 2), [1, 4]);
        msh.faceNums( 1, faceChildren(:, 6)) = 6;         
    end 
    % 12 interior faces
    tempIDs = msh.nFaces + 1 : msh.nFaces + 4;
    faceIDs1(2, :, :) = reshape(tempIDs, [1, 2, 2]);
    msh.faceElems(1, tempIDs) = reshape(elemChildren(1, :, :), [1, 4]);
    msh.faceElems(2, tempIDs) = reshape(elemChildren(2, :, :), [1, 4]);
    msh.faceNums( 1, tempIDs) = 2;
    msh.faceNums( 2, tempIDs) = 1;
    msh.faceNormalx(tempIDs)  = 1;
    msh.faceSize(tempIDs)     = h(2) * h(3) / 4;
    msh.faceJac(tempIDs)      = h(2) * h(3) / 16;
    msh.nFaces = msh.nFaces + 4; 
    
    tempIDs = msh.nFaces + 1 : msh.nFaces + 4;
    faceIDs2(:, 2, :) = reshape(msh.nFaces + [1, 3, 2, 4], [2, 1, 2]); 
    msh.faceElems(1, tempIDs) = elemChildren([1, 5, 2, 6]);
    msh.faceElems(2, tempIDs) = elemChildren([3, 7, 4, 8]);
    msh.faceNums( 1, tempIDs) = 4;
    msh.faceNums( 2, tempIDs) = 3;
    msh.faceNormaly(tempIDs)  = 1;
    msh.faceSize(tempIDs)     = h(1) * h(3) / 4;
    msh.faceJac(tempIDs)      = h(1) * h(3) / 16;
    msh.nFaces = msh.nFaces + 4;   
    
    tempIDs = msh.nFaces + 1 : msh.nFaces + 4;
    faceIDs3(:, :, 2) = reshape(tempIDs, [2, 2]);
    msh.faceElems(1, tempIDs) = reshape(elemChildren(:, :, 1), [1, 4]);
    msh.faceElems(2, tempIDs) = reshape(elemChildren(:, :, 2), [1, 4]);
    msh.faceNums( 1, tempIDs) = 6;
    msh.faceNums( 2, tempIDs) = 5;
    msh.faceNormalz(tempIDs)  = 1;
    msh.faceSize(tempIDs)     = h(1) * h(2) / 4;
    msh.faceJac(tempIDs)      = h(1) * h(2) / 16;
    msh.nFaces = msh.nFaces + 4;    
    
    % element-face connections
    msh.elemFaces(1, elemChildren) = reshape(faceIDs1(1 : 2, :, :), [1, 8]);
    msh.elemFaces(2, elemChildren) = reshape(faceIDs1(2 : 3, :, :), [1, 8]);
    msh.elemFaces(3, elemChildren) = reshape(faceIDs2(:, 1 : 2, :), [1, 8]);
    msh.elemFaces(4, elemChildren) = reshape(faceIDs2(:, 2 : 3, :), [1, 8]);
    msh.elemFaces(5, elemChildren) = reshape(faceIDs3(:, :, 1 : 2), [1, 8]);
    msh.elemFaces(6, elemChildren) = reshape(faceIDs3(:, :, 2 : 3), [1, 8]);
    msh.nElems = msh.nElems + 8;                 
    
    % Project the solution in the coarse element to the eight refined
    % elements
    UE = UT(:, nv * (ID - 1) + 1 : nv * ID);
    UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(8)) = [c2f{1} * UE, c2f{2} * UE, c2f{3} * UE, c2f{4} * UE, c2f{5} * UE, c2f{6} * UE, c2f{7} * UE, c2f{8} * UE];    
end
msh.faceNormalx(:, msh.nFaces + 1 : end)  = []; 
msh.faceNormaly(:, msh.nFaces + 1 : end)  = []; 
msh.faceNormalz(:, msh.nFaces + 1 : end)  = [];
msh.faceSize(msh.nFaces + 1 : end)        = []; 
msh.faceType(msh.nFaces + 1 : end)        = []; 
msh.faceElems(:, msh.nFaces + 1 : end)    = []; 
msh.faceNums(:, msh.nFaces + 1 : end)     = []; 
msh.faceLevel(msh.nFaces + 1 : end)       = []; 
msh.faceParent(msh.nFaces + 1 : end)      = []; 
msh.faceChildren(:, msh.nFaces + 1 : end) = []; 
msh.faceJac(:, msh.nFaces + 1 : end)      = [];

%**************************************************************************
% Then coarsen the elements. Note that we may delete some elements to
% coarsen to ensure the level difference between neighboring elements is no 
% more than 1
%**************************************************************************
% Delete the candidate elements which don't have a parent or should refine
ceIDs = candCeIDs(msh.elemParent(candCeIDs) ~= 0 & ~shouldRefine(candCeIDs));
ceIDs = sort(ceIDs);
nce = length(ceIDs);
shouldCoarsen = false(1, msh.nElems);
shouldCoarsen(ceIDs) = true;

% Make sure all the eight siblings need to be coarsened and store the IDs of
% their parents
npe = 0;
peIDs = zeros(1, floor(nce / 8));
i = 1;
while i <= nce
    ID = ceIDs(i);
    parentID = msh.elemParent(ID);
    % Make sure all the eight siblings need to be coarsened
    siblingIDs = msh.elemChildren(:, parentID);  
    if any(~shouldCoarsen(siblingIDs))
        shouldCoarsen(ID) = false;
        i = i + 1;
        continue;
    end
    % Store the parent ID
    npe = npe + 1;
    peIDs(npe) = parentID;
    
    i = i + 8;
end
peIDs(npe + 1 : end) = [];

% We should first coarsen the elements of level maxLevel, then level 
% maxLevel - 1, then ..., finally level 1
peLevels = msh.elemLevel(peIDs);
speIDs = zeros(1, npe);
s = 0;
for i = msh.maxLevel - 1 : -1 : 0
    temp = peIDs(peLevels == i);
    l = length(temp);
    speIDs(s + 1 : s + l) = temp;
    s = s + l;
end
if (s ~= npe)
    error('Wrong size after sorting')
end
peIDs = speIDs;

% Update face-element connections, or delete subelements and subfaces and 
% store their IDs 
deIDs = zeros(1, 8 * npe);
dfIDs = zeros(1, 36 * npe);
nde = 0;
ndf = 0;
for i = 1 : npe
    ID = peIDs(i);
    elemFaces    = msh.elemFaces(:, ID);
    elemChildren = msh.elemChildren(:, ID); 
    faceTypes    = msh.faceType(elemFaces);
    faceElems    = msh.faceElems(:, elemFaces);
    faceChildren = msh.faceChildren(:, elemFaces);
    
    % Make sure the level difference of adjacent elements is less than 1
    if any(msh.faceChildren(1, faceChildren) ~= 0)
        shouldCoarsen(elemChildren) = false;
        continue;  
    end
    
    % Project the solution in the eight refined elements into the coarse
    % element
    UT(:, nv * (ID - 1) + 1 : nv * ID) = f2c{1} * UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) ... 
                                       + f2c{2} * UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)) ...
                                       + f2c{3} * UT(:, nv * (elemChildren(3) - 1) + 1 : nv * elemChildren(3)) ...
                                       + f2c{4} * UT(:, nv * (elemChildren(4) - 1) + 1 : nv * elemChildren(4)) ...
                                       + f2c{5} * UT(:, nv * (elemChildren(5) - 1) + 1 : nv * elemChildren(5)) ...
                                       + f2c{6} * UT(:, nv * (elemChildren(6) - 1) + 1 : nv * elemChildren(6)) ...
                                       + f2c{7} * UT(:, nv * (elemChildren(7) - 1) + 1 : nv * elemChildren(7)) ...
                                       + f2c{8} * UT(:, nv * (elemChildren(8) - 1) + 1 : nv * elemChildren(8));
    
    % Update face-element connections or delete subfaces
    % first face
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        ID1 = faceElems(1, 1);
        if (msh.elemChildren(1, ID1) ~= 0)
            msh.faceElems(2, faceChildren(:, 1)) = ID;
            msh.faceNums( 2, faceChildren(:, 1)) = 7 : 10;            
        else
            msh.faceChildren(:, elemFaces(1)) = 0;
            dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 1)';
            ndf = ndf + 4;            
        end
    else
        msh.faceChildren(:, elemFaces(1)) = 0;    
        dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 1)';
        ndf = ndf + 4;         
    end
    % second face
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        ID2 = faceElems(2, 2);
        if (msh.elemChildren(1, ID2) ~= 0)
            msh.faceElems(1, faceChildren(:, 2)) = ID;
            msh.faceNums( 1, faceChildren(:, 2)) = 11 : 14;             
        else
            msh.faceChildren(:, elemFaces(2)) = 0;
            dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 2)';
            ndf = ndf + 4;            
        end        
    else
        msh.faceChildren(:, elemFaces(2)) = 0;
        dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 2)';
        ndf = ndf + 4;          
    end
    % third face
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        ID3 = faceElems(1, 3);
        if (msh.elemChildren(1, ID3) ~= 0)
            msh.faceElems(2, faceChildren(:, 3)) = ID;
            msh.faceNums( 2, faceChildren(:, 3)) = 15 : 18;             
        else
            msh.faceChildren(:, elemFaces(3)) = 0;
            dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 3)';
            ndf = ndf + 4;  
        end        
    else
        msh.faceChildren(:, elemFaces(3)) = 0;
        dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 3)';
        ndf = ndf + 4;         
    end
    % fourth face
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        ID4 = faceElems(2, 4);
        if (msh.elemChildren(1, ID4) ~= 0)
            msh.faceElems(1, faceChildren(:, 4)) = ID;
            msh.faceNums( 1, faceChildren(:, 4)) = 19 : 22;             
        else
            msh.faceChildren(:, elemFaces(4)) = 0;
            dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 4)';
            ndf = ndf + 4; 
        end        
    else
        msh.faceChildren(:, elemFaces(4)) = 0;
        dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 4)';
        ndf = ndf + 4; 
    end
    % fifth face
    if (faceTypes(5) == 0 || faceTypes(5) == 1)
        ID5 = faceElems(1, 5);
        if (msh.elemChildren(1, ID5) ~= 0)
            msh.faceElems(2, faceChildren(:, 5)) = ID;
            msh.faceNums( 2, faceChildren(:, 5)) = 23 : 26;             
        else
            msh.faceChildren(:, elemFaces(5)) = 0;
            dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 5)';
            ndf = ndf + 4;             
        end        
    else
        msh.faceChildren(:, elemFaces(5)) = 0;
        dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 5)';
        ndf = ndf + 4;         
    end
    % sixth face
    if (faceTypes(6) == 0 || faceTypes(6) == 1)
        ID6 = faceElems(2, 6);
        if (msh.elemChildren(1, ID6) ~= 0)
            msh.faceElems(1, faceChildren(:, 6)) = ID;
            msh.faceNums( 1, faceChildren(:, 6)) = 27 : 30;             
        else
            msh.faceChildren(:, elemFaces(6)) = 0;
            dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 6)';
            ndf = ndf + 4;             
        end        
    else
        msh.faceChildren(:, elemFaces(6)) = 0;
        dfIDs(ndf + 1 : ndf + 4) = faceChildren(:, 6)';
        ndf = ndf + 4;         
    end
    % 12 interior faces
    childFaceID = msh.elemFaces(2, elemChildren(1));
    dfIDs(ndf + 1 : ndf + 12) = childFaceID : childFaceID + 11;
    ndf = ndf + 12;
    
    % Delete elements
    msh.elemChildren(:, ID) = 0;
    deIDs(nde + 1 : nde + 8) = elemChildren;
    nde = nde + 8;
end
deIDs(nde + 1 : end) = [];
dfIDs(ndf + 1 : end) = [];

% Re-index the remaining elements and faces
remEIDs = 1 : msh.nElems;
remFIDs = 1 : msh.nFaces;
remEIDs(deIDs) = []; % remaining element old IDs
remFIDs(dfIDs) = []; % remaining face old IDs

neIDs = zeros(1, msh.nElems); 
nfIDs = zeros(1, msh.nFaces); 
neIDs(remEIDs) = 1 : length(remEIDs); % new element IDs
nfIDs(remFIDs) = 1 : length(remFIDs); % new face IDs

% Refresh element-face connections, parent and child IDs of elements
msh.elemFaces(:, remEIDs) = nfIDs(msh.elemFaces(:, remEIDs));

remEIDs_hasParent = remEIDs(msh.elemParent(remEIDs) ~= 0);
msh.elemParent(remEIDs_hasParent) = neIDs(msh.elemParent(remEIDs_hasParent));

remEIDs_hasChild = remEIDs(msh.elemChildren(1, remEIDs) ~= 0);
msh.elemChildren(:, remEIDs_hasChild) = neIDs(msh.elemChildren(:, remEIDs_hasChild));

% Refresh face-element connections, parent and child IDs of faces
msh.faceElems(1, remFIDs) = neIDs(msh.faceElems(1, remFIDs));
remFIDs_int = remFIDs(msh.faceType(remFIDs) == 0 | msh.faceType(remFIDs) == 1);
msh.faceElems(2, remFIDs_int) = neIDs(msh.faceElems(2, remFIDs_int));

remFIDs_hasParent = remFIDs(msh.faceParent(remFIDs) ~= 0);
msh.faceParent(remFIDs_hasParent) = nfIDs(msh.faceParent(remFIDs_hasParent));

remFIDs_hasChild = remFIDs(msh.faceChildren(1, remFIDs) ~= 0);
msh.faceChildren(:, remFIDs_hasChild) = nfIDs(msh.faceChildren(:, remFIDs_hasChild));

% Finally remove the deleted elements and faces from the mesh
msh.elemCenter(:, deIDs)   = [];
msh.elemLength(:, deIDs)   = [];
msh.elemSize(deIDs)        = [];
msh.elemFaces(:, deIDs)    = [];
msh.elemLevel(deIDs)       = [];
msh.elemParent(deIDs)      = [];
msh.elemChildren(:, deIDs) = [];
msh.elemJac(:, deIDs)      = [];
msh.elemJxix(:, deIDs)     = [];
msh.elemJetay(:, deIDs)    = [];
msh.elemJzetaz(:, deIDs)   = [];
msh.faceNormalx(:, dfIDs)  = [];
msh.faceNormaly(:, dfIDs)  = [];
msh.faceNormalz(:, dfIDs)  = [];
msh.faceSize(dfIDs)        = [];
msh.faceType(dfIDs)        = [];
msh.faceElems(:, dfIDs)    = [];
msh.faceNums(:, dfIDs)     = [];
msh.faceLevel(dfIDs)       = []; 
msh.faceParent(dfIDs)      = []; 
msh.faceChildren(:, dfIDs) = []; 
msh.faceJac(:, dfIDs)      = [];
msh.nElems = length(remEIDs);
msh.nFaces = length(remFIDs);

%**************************************************************************
% Finally update the IDs of leaf elements and faces and get the new
% solution
%**************************************************************************
% Update the IDs of leaf elements
msh.LElems  = find(msh.elemChildren(1, :) == 0);
msh.nLElems = length(msh.LElems);

msh.elemLID             = zeros(1, msh.nElems);
msh.elemLID(msh.LElems) = 1 : msh.nLElems;

% Update the IDs of leaf faces
msh.LFaces  = find(msh.faceChildren(1, :) == 0);
msh.nLFaces = length(msh.LFaces);

msh.faceLID             = zeros(1, msh.nFaces);
msh.faceLID(msh.LFaces) = 1 : msh.nLFaces;

% Update the IDs of internal and various boundary faces
isInternal     = msh.faceType(msh.LFaces) == 0;
msh.intLFaces  = msh.LFaces(isInternal);
msh.nIntLFaces = length(msh.intLFaces);
bndLFaces      = msh.LFaces(~isInternal);
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = bndLFaces(msh.faceType(bndLFaces) == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end

% Update the solution in the new leaf elements
UT(:, (-nv + 1 : 0)' + nv * deIDs) = [];
U = UT(:, (-nv + 1 : 0)' + nv * msh.LElems);
if (layout == 1)
    U = U(:, (1 : nv : end)' + (0 : nv - 1));
end

end








