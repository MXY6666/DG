% Suppose the level difference between neighboring elements is no more than
% 1 when refining and coarsening
% Copyright (C) Xiangyi Meng.

function [msh, U] = refineAndCoarsenRectMesh(msh, candReIDs, candCeIDs, nv, U, c2f, f2c, layout)

if (nargin < 1)
    error('Not enough arguments')
end

if (msh.type ~= 202) 
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
    U = zeros(3, msh.nLElems);
end
if (size(U, 2) ~= msh.nLElems * nv)
    error('The given mesh and the solution coefficients does not coinside')
end

if (nargin < 6) || isempty(c2f)
    c2f = getC2FMatrix_rect(1);
end

if (nargin < 7) || isempty(f2c)
    f2c = getF2CMatrix_rect(1);
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
    
    % Determine if refine the left adjacent element
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        IDL = faceElems(1, 1);
        if (msh.elemLevel(IDL) < level) && ~shouldRefine(IDL)
            shouldRefine(IDL) = true;
            nre = nre + 1;
            reIDs(nre) = IDL;
        end         
    end
    % Determine if refine the right adjacent element
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        IDR = faceElems(2, 2);
        if (msh.elemLevel(IDR) < level) && ~shouldRefine(IDR)
            shouldRefine(IDR) = true;
            nre = nre + 1;
            reIDs(nre) = IDR;
        end        
    end
    % Determine if refine the bottom adjacent element
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        IDD = faceElems(1, 3);
        if (msh.elemLevel(IDD) < level) && ~shouldRefine(IDD)
            shouldRefine(IDD) = true;
            nre = nre + 1;
            reIDs(nre) = IDD;
        end        
    end    
    % Determine if refine the top adjacent element
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        IDU = faceElems(2, 4);
        if (msh.elemLevel(IDU) < level) && ~shouldRefine(IDU)
            shouldRefine(IDU) = true;
            nre = nre + 1;
            reIDs(nre) = IDU;
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
msh.elemCenter(:, msh.nElems + 1 : msh.nElems + 4 * nre)   = 0;
msh.elemLength(:, msh.nElems + 1 : msh.nElems + 4 * nre)   = 0;
msh.elemSize(msh.nElems + 1 : msh.nElems + 4 * nre)        = 0;
msh.elemFaces(:, msh.nElems + 1 : msh.nElems + 4 * nre)    = 0;
msh.elemLevel(msh.nElems + 1 : msh.nElems + 4 * nre)       = 0;
msh.elemParent(msh.nElems + 1 : msh.nElems + 4 * nre)      = 0;
msh.elemChildren(:, msh.nElems + 1 : msh.nElems + 4 * nre) = 0;
msh.elemJac(:, msh.nElems + 1 : msh.nElems + 4 * nre)      = 0;
msh.elemJxix(:, msh.nElems + 1 : msh.nElems + 4 * nre)     = 0;
msh.elemJetay(:, msh.nElems + 1 : msh.nElems + 4 * nre)    = 0;
% preallocation for faces
msh.faceNormalx(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)  = 0; 
msh.faceNormaly(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)  = 0; 
msh.faceSize(msh.nFaces + 1 : msh.nFaces + 12 * nre)        = 0; 
msh.faceType(msh.nFaces + 1 : msh.nFaces + 12 * nre)        = 0; 
msh.faceElems(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)    = 0; 
msh.faceNums(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)     = 0; 
msh.faceLevel(msh.nFaces + 1 : msh.nFaces + 12 * nre)       = 0; 
msh.faceParent(msh.nFaces + 1 : msh.nFaces + 12 * nre)      = 0; 
msh.faceChildren(:, msh.nFaces + 1 : msh.nFaces + 12 * nre) = 0; 
msh.faceJac(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)      = 0; 
% preallocation for solution
UT = zeros(size(U, 1), (msh.nElems + 4 * nre) * nv);
if (layout == 1)
    U = U(:, (0 : nv - 1)' * msh.nLElems + (1 : msh.nLElems));
end
UT(:, (-nv + 1 : 0)' + nv * msh.LElems) = U;


subfIDs = zeros(2, 6);
for i = 1 : nre
    % current loop element and its faces
    ID = reIDs(i);
    ct = msh.elemCenter(:, ID);
    h  = msh.elemLength(:, ID);
    si = msh.elemSize(ID);
    elemFaces    = msh.elemFaces(:, ID);
    faceNormalxs = msh.faceNormalx(elemFaces);
    faceNormalys = msh.faceNormaly(elemFaces);
    faceTypes    = msh.faceType(elemFaces);
    faceElems    = msh.faceElems(:, elemFaces);
    faceLevels   = msh.faceLevel(elemFaces);
    faceChildren = msh.faceChildren(:, elemFaces);

    % creation of subelements
    elemChildren = msh.nElems + 1 : msh.nElems + 4;
    msh.elemChildren(:, ID)         = elemChildren';    
    msh.elemCenter(:, elemChildren) = [ct(1) - h(1) / 4, ct(1) + h(1) / 4, ct(1) - h(1) / 4, ct(1) + h(1) / 4; ...
                                       ct(2) - h(2) / 4, ct(2) - h(2) / 4, ct(2) + h(2) / 4, ct(2) + h(2) / 4];
    msh.elemLength(:, elemChildren) = repmat(h / 2, [1, 4]);
    msh.elemSize(elemChildren)      = si / 4; 
    msh.elemLevel(elemChildren)     = msh.elemLevel(ID) + 1;
    msh.elemParent(elemChildren)    = ID;
    msh.elemJac(elemChildren)       = si / 16;
    msh.elemJxix(elemChildren)      = h(2) / 4;
    msh.elemJetay(elemChildren)     = h(1) / 4;

    % creation of subfaces
    % left face
    if (faceChildren(1, 1) == 0)
        msh.faceChildren(:, elemFaces(1)) = [msh.nFaces + 1; msh.nFaces + 2];

        subfIDs(:, 1) = [msh.nFaces + 1; msh.nFaces + 2];
        if (faceTypes(1) == 0 || faceTypes(1) == 1)
            msh.faceElems(1, subfIDs(:, 1)) = faceElems(1, 1);
            msh.faceElems(2, subfIDs(:, 1)) = elemChildren([1, 3]);
            msh.faceNums( 1, subfIDs(:, 1)) = [7, 8];
            msh.faceNums( 2, subfIDs(:, 1)) = 1;
        else
            msh.faceElems(1, subfIDs(:, 1)) = elemChildren([1, 3]);
            msh.faceNums( 1, subfIDs(:, 1)) = 1;
        end   
        msh.faceNormalx(subfIDs(:, 1)) = faceNormalxs(1);
        msh.faceSize(subfIDs(:, 1))    = h(2) / 2;
        msh.faceType(subfIDs(:, 1))    = faceTypes(1);
        msh.faceLevel(subfIDs(:, 1))   = faceLevels(1) + 1;
        msh.faceParent(subfIDs(:, 1))  = elemFaces(1);
        msh.faceJac(subfIDs(:, 1))     = h(2) / 4;
        msh.nFaces = msh.nFaces + 2;        
    else
        subfIDs(:, 1) = faceChildren(:, 1); 
        % It must be an internal or periodic face in this case
        msh.faceElems(2, subfIDs(:, 1)) = elemChildren([1, 3]);
        msh.faceNums( 2, subfIDs(:, 1)) = 1;
    end
    % right face
    if (faceChildren(1, 2) == 0)
        msh.faceChildren(:, elemFaces(2)) = [msh.nFaces + 1; msh.nFaces + 2];

        subfIDs(:, 2) = [msh.nFaces + 1; msh.nFaces + 2];
        msh.faceElems(1, subfIDs(:, 2)) = elemChildren([2, 4]);
        msh.faceNums( 1, subfIDs(:, 2)) = 2;
        if (faceTypes(2) == 0 || faceTypes(2) == 1)
            msh.faceElems(2, subfIDs(:, 2)) = faceElems(2, 2);
            msh.faceNums( 2, subfIDs(:, 2)) = [5, 6];
        end
        msh.faceNormalx(subfIDs(:, 2)) = faceNormalxs(2);
        msh.faceSize(subfIDs(:, 2))    = h(2) / 2;
        msh.faceType(subfIDs(:, 2))    = faceTypes(2);
        msh.faceLevel(subfIDs(:, 2))   = faceLevels(2) + 1;
        msh.faceParent(subfIDs(:, 2))  = elemFaces(2);
        msh.faceJac(subfIDs(:, 2))     = h(2) / 4;
        msh.nFaces = msh.nFaces + 2;        
    else
        subfIDs(:, 2) = faceChildren(:, 2); 
        % It must be an internal or periodic face in this case
        msh.faceElems(1, subfIDs(:, 2)) = elemChildren([2, 4]);
        msh.faceNums( 1, subfIDs(:, 2)) = 2;
    end   
    % bottom face
    if (faceChildren(1, 3) == 0)
        msh.faceChildren(:, elemFaces(3)) = [msh.nFaces + 1; msh.nFaces + 2];

        subfIDs(:, 3) = [msh.nFaces + 1; msh.nFaces + 2];
        if (faceTypes(3) == 0 || faceTypes(3) == 1)
            msh.faceElems(1, subfIDs(:, 3)) = faceElems(1, 3);
            msh.faceElems(2, subfIDs(:, 3)) = elemChildren([1, 2]);
            msh.faceNums( 1, subfIDs(:, 3)) = [11, 12];
            msh.faceNums( 2, subfIDs(:, 3)) = 3;
        else
            msh.faceElems(1, subfIDs(:, 3)) = elemChildren([1, 2]);
            msh.faceNums( 1, subfIDs(:, 3)) = 3;
        end     
        msh.faceNormaly(subfIDs(:, 3)) = faceNormalys(3);
        msh.faceSize(subfIDs(:, 3))    = h(1) / 2;
        msh.faceType(subfIDs(:, 3))    = faceTypes(3);
        msh.faceLevel(subfIDs(:, 3))   = faceLevels(3) + 1;
        msh.faceParent(subfIDs(:, 3))  = elemFaces(3);
        msh.faceJac(subfIDs(:, 3))     = h(1) / 4;
        msh.nFaces = msh.nFaces + 2;
    else
        subfIDs(:, 3) = faceChildren(:, 3); 
        % It must be an internal or periodic face in this case
        msh.faceElems(2, subfIDs(:, 3)) = elemChildren([1, 2]);
        msh.faceNums( 2, subfIDs(:, 3)) = 3;    
    end    
    % top face
    if (faceChildren(1, 4) == 0)
        msh.faceChildren(:, elemFaces(4)) = [msh.nFaces + 1; msh.nFaces + 2];

        subfIDs(:, 4) = [msh.nFaces + 1; msh.nFaces + 2]; 
        msh.faceElems(1, subfIDs(:, 4)) = elemChildren([3, 4]);
        msh.faceNums( 1, subfIDs(:, 4)) = 4;
        if (faceTypes(4) == 0 || faceTypes(4) == 1)
            msh.faceElems(2, subfIDs(:, 4)) = faceElems(2, 4);  
            msh.faceNums( 2, subfIDs(:, 4)) = [9, 10];
        end
        msh.faceNormaly(subfIDs(:, 4)) = faceNormalys(4);
        msh.faceSize(subfIDs(:, 4))    = h(1) / 2;
        msh.faceType(subfIDs(:, 4))    = faceTypes(4);
        msh.faceLevel(subfIDs(:, 4))   = faceLevels(4) + 1;
        msh.faceParent(subfIDs(:, 4))  = elemFaces(4);
        msh.faceJac(subfIDs(:, 4))     = h(1) / 4;
        msh.nFaces = msh.nFaces + 2;        
    else
        subfIDs(:, 4) = faceChildren(:, 4);  
        % It must be an internal or periodic face in this case
        msh.faceElems(1, subfIDs(:, 4)) = elemChildren([3, 4]);
        msh.faceNums( 1, subfIDs(:, 4)) = 4;
    end   
    % four interior faces
    subfIDs(:, 5) = [msh.nFaces + 1; msh.nFaces + 2]; 
    msh.faceElems(1, subfIDs(:, 5)) = elemChildren([1, 3]);
    msh.faceElems(2, subfIDs(:, 5)) = elemChildren([2, 4]);
    msh.faceNums( 1, subfIDs(:, 5)) = 2;
    msh.faceNums( 2, subfIDs(:, 5)) = 1;
    msh.faceNormalx(subfIDs(:, 5))  = 1;
    msh.faceSize(subfIDs(:, 5))     = h(2) / 2;
    msh.faceJac(subfIDs(:, 5))      = h(2) / 4;
    subfIDs(:, 6) = [msh.nFaces + 3; msh.nFaces + 4]; 
    msh.faceElems(1, subfIDs(:, 6)) = elemChildren([1, 2]);
    msh.faceElems(2, subfIDs(:, 6)) = elemChildren([3, 4]);
    msh.faceNums( 1, subfIDs(:, 6)) = 4;
    msh.faceNums( 2, subfIDs(:, 6)) = 3;
    msh.faceNormaly(subfIDs(:, 6))  = 1;
    msh.faceSize(subfIDs(:, 6))     = h(1) / 2; 
    msh.faceJac(subfIDs(:, 6))      = h(1) / 4; 
    msh.nFaces = msh.nFaces + 4;
    
    % element-face connections
    msh.elemFaces(:, elemChildren) = [subfIDs(1, 1), subfIDs(1, 5), subfIDs(2, 1), subfIDs(2, 5); ...
                                      subfIDs(1, 5), subfIDs(1, 2), subfIDs(2, 5), subfIDs(2, 2); ...
                                      subfIDs(1, 3), subfIDs(2, 3), subfIDs(1, 6), subfIDs(2, 6); ...
                                      subfIDs(1, 6), subfIDs(2, 6), subfIDs(1, 4), subfIDs(2, 4)];
    msh.nElems = msh.nElems + 4;

    % Project the solution in the coarse element to the four refined
    % elements
    UE = UT(:, nv * (ID - 1) + 1 : nv * ID);
    UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(4)) = [c2f{1} * UE, c2f{2} * UE, c2f{3} * UE, c2f{4} * UE];
end    
msh.faceNormalx(:, msh.nFaces + 1 : end)  = []; 
msh.faceNormaly(:, msh.nFaces + 1 : end)  = []; 
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

% Make sure all the four siblings need to be coarsened and store the IDs of
% their parents
npe = 0;
peIDs = zeros(1, floor(nce / 4));
i = 1;
while i <= nce
    ID = ceIDs(i);
    parentID = msh.elemParent(ID);
    % Make sure all the four siblings need to be coarsened
    siblingIDs = msh.elemChildren(:, parentID); 
    if any(~shouldCoarsen(siblingIDs))
        shouldCoarsen(ID) = false;
        i = i + 1;
        continue;
    end
    % Store the parent ID
    npe = npe + 1;
    peIDs(npe) = parentID;
    
    i = i + 4;
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
deIDs = zeros(1, 4 * npe);
dfIDs = zeros(1, 12 * npe);
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
    
    % Project the solution in the four refined elements into the coarse
    % element
    UT(:, nv * (ID - 1) + 1 : nv * ID) = f2c{1} * UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) ...
                                       + f2c{2} * UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)) ...
                                       + f2c{3} * UT(:, nv * (elemChildren(3) - 1) + 1 : nv * elemChildren(3)) ...
                                       + f2c{4} * UT(:, nv * (elemChildren(4) - 1) + 1 : nv * elemChildren(4));
                        
    % Update face-element connections or delete subfaces
    % left face
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        IDL = faceElems(1, 1);
        if (msh.elemChildren(1, IDL) ~= 0)
            msh.faceElems(2, faceChildren(:, 1)) = ID;
            msh.faceNums( 2, faceChildren(:, 1)) = [5, 6];
        else
            msh.faceChildren(:, elemFaces(1)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 1)';
            ndf = ndf + 2;           
        end
    else
        msh.faceChildren(:, elemFaces(1)) = 0;       
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 1)';
        ndf = ndf + 2;       
    end   
    % right face
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        IDR = faceElems(2, 2);
        if (msh.elemChildren(1, IDR) ~= 0)
            msh.faceElems(1, faceChildren(:, 2)) = ID;
            msh.faceNums( 1, faceChildren(:, 2)) = [7, 8];
        else
            msh.faceChildren(:, elemFaces(2)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 2)';
            ndf = ndf + 2;             
        end
    else
        msh.faceChildren(:, elemFaces(2)) = 0;       
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 2)';
        ndf = ndf + 2;
    end  
    % bottom face
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        IDD = faceElems(1, 3);
        if (msh.elemChildren(1, IDD) ~= 0)
            msh.faceElems(2, faceChildren(:, 3)) = ID;
            msh.faceNums( 2, faceChildren(:, 3)) = [9, 10];
        else
            msh.faceChildren(:, elemFaces(3)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 3)';
            ndf = ndf + 2;
        end
    else        
        msh.faceChildren(:, elemFaces(3)) = 0;
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 3)';
        ndf = ndf + 2;
    end    
    % top face
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        IDU = faceElems(2, 4);
        if (msh.elemChildren(1, IDU) ~= 0)
            msh.faceElems(1, faceChildren(:, 4)) = ID;
            msh.faceNums( 1, faceChildren(:, 4)) = [11, 12];
        else          
            msh.faceChildren(:, elemFaces(4)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 4)';
            ndf = ndf + 2;           
        end
    else
        msh.faceChildren(:, elemFaces(4)) = 0;        
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 4)';
        ndf = ndf + 2;        
    end  
    % four interior faces
    childFaceID = msh.elemFaces(2, elemChildren(1));
    dfIDs(ndf + 1 : ndf + 4) = childFaceID : childFaceID + 3;
    ndf = ndf + 4;
    
    % Delete elements
    msh.elemChildren(:, ID) = 0;
    deIDs(nde + 1 : nde + 4) = elemChildren;
    nde = nde + 4;
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
msh.faceNormalx(:, dfIDs)  = [];
msh.faceNormaly(:, dfIDs)  = [];
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




