% Suppose the level difference between neighboring elements is no more than
% 1 when refining and coarsening
% Copyright (C) Xiangyi Meng.

function [msh, U] = refineAndCoarsenLineMesh(msh, candReIDs, candCeIDs, nv, U, c2f, f2c, layout)

if (nargin < 1)
    error('Not enough arguments')
end

if (msh.type ~= 102) 
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
    U = zeros(2, msh.nLElems);
end
if (size(U, 2) ~= msh.nLElems * nv)
    error('The given mesh and the solution coefficients does not coinside')
end

if (nargin < 6) || isempty(c2f)
    c2f = getC2FMatrix_line(1);
end

if (nargin < 7) || isempty(f2c)
    f2c = getF2CMatrix_line(1);
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
    i = i + 1;
end
reIDs(nre + 1 : end) = [];

% We refine the elements by adding new elements and faces
% preallocation for elements
msh.elemCenter(:, msh.nElems + 1 : msh.nElems + 2 * nre)   = 0;
msh.elemLength(:, msh.nElems + 1 : msh.nElems + 2 * nre)   = 0;
msh.elemSize(msh.nElems + 1 : msh.nElems + 2 * nre)        = 0;
msh.elemFaces(:, msh.nElems + 1 : msh.nElems + 2 * nre)    = 0;
msh.elemLevel(msh.nElems + 1 : msh.nElems + 2 * nre)       = 0;
msh.elemParent(msh.nElems + 1 : msh.nElems + 2 * nre)      = 0;
msh.elemChildren(:, msh.nElems + 1 : msh.nElems + 2 * nre) = 0;
msh.elemJac(:, msh.nElems + 1 : msh.nElems + 2 * nre)      = 0;
% preallocation for faces
msh.faceNormalx(:, msh.nFaces + 1 : msh.nFaces + nre)  = 1;  
msh.faceType(msh.nFaces + 1 : msh.nFaces + nre)        = 0; 
msh.faceElems(:, msh.nFaces + 1 : msh.nFaces + nre)    = 0; 
msh.faceNums(:, msh.nFaces + 1 : msh.nFaces + nre)     = repmat([2; 1], [1, nre]);  
% preallocation for solution
UT = zeros(size(U, 1), (msh.nElems + 2 * nre) * nv);
if (layout == 1)
    U = U(:, (0 : nv - 1)' * msh.nLElems + (1 : msh.nLElems));
end
UT(:, (-nv + 1 : 0)' + nv * msh.LElems) = U;

for i = 1 : nre
    % current loop element and its faces
    ID = reIDs(i);
    ct = msh.elemCenter(:, ID);
    h  = msh.elemLength(:, ID);
    elemFaces = msh.elemFaces(:, ID);
   
    % creation of subelements
    elemChildren = msh.nElems + 1 : msh.nElems + 2;
    msh.elemChildren(:, ID)         = elemChildren'; 
    msh.elemCenter(:, elemChildren) = [ct - h / 4, ct + h / 4];
    msh.elemLength(:, elemChildren) = h / 2;
    msh.elemSize(elemChildren)      = h / 2;
    msh.elemLevel(elemChildren)     = msh.elemLevel(ID) + 1;
    msh.elemParent(elemChildren)    = ID;
    msh.elemJac(elemChildren)       = h / 4;
    
    % Update the information of existing faces
    msh.faceElems(2, elemFaces(1)) = elemChildren(1);
    msh.faceElems(1, elemFaces(2)) = elemChildren(2);
    
    % creation of new faces
    faceID = msh.nFaces + 1;
    msh.faceElems(:, faceID) = elemChildren'; 
    msh.nFaces = msh.nFaces + 1;

    % element-face connections
    msh.elemFaces(:, elemChildren) = [elemFaces(1), faceID; faceID, elemFaces(2)]; 
    msh.nElems = msh.nElems + 2;
    
    % Project the solution in the coarse element to the two refined
    % elements
    UE = UT(:, nv * (ID - 1) + 1 : nv * ID);
    UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(2)) = [c2f{1} * UE, c2f{2} * UE];    
end

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

% Make sure both the siblings need to be coarsened and store the IDs of
% their parents
npe = 0;
peIDs = zeros(1, floor(nce / 2));
i = 1;
while i <= nce
    ID = ceIDs(i);
    parentID = msh.elemParent(ID);
    % Make sure both the siblings need to be coarsened
    siblingIDs = msh.elemChildren(:, parentID); 
    if any(~shouldCoarsen(siblingIDs))
        shouldCoarsen(ID) = false;
        i = i + 1;
        continue;
    end
    % Store the parent ID
    npe = npe + 1;
    peIDs(npe) = parentID;
    
    i = i + 2;
end
peIDs(npe + 1 : end) = [];

% Update face-element connections, or delete subelements and faces and 
% store their IDs 
deIDs = zeros(1, 2 * npe);
dfIDs = zeros(1, npe);
nde = 0;
ndf = 0;
for i = 1 : npe
    ID = peIDs(i);
    elemFaces    = msh.elemFaces(:, ID);
    elemLevel    = msh.elemLevel(ID);
    elemChildren = msh.elemChildren(:, ID); 
    faceElems    = msh.faceElems(:, elemFaces);
    eleLevelL    = msh.elemLevel(faceElems(1, 1));
    eleLevelR    = msh.elemLevel(faceElems(2, 2));
    
    % Make sure the level difference of adjacent elements is less than 1
    if (eleLevelL > elemLevel + 1) || (eleLevelR > elemLevel + 1)
        shouldCoarsen(elemChildren) = false;
        continue;        
    end
    
    % Project the solution in the two refined elements into the coarse
    % element
    UT(:, nv * (ID - 1) + 1 : nv * ID) = f2c{1} * UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) ...
                                       + f2c{2} * UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)); 
                                   
    % Update face-element connections and delete faces
    msh.faceElems(2, elemFaces(1)) = ID;
    msh.faceElems(1, elemFaces(2)) = ID;
    dfIDs(ndf + 1) = msh.elemFaces(2, elemChildren(1));
    ndf = ndf + 1;  
 
    % Delete elements
    msh.elemChildren(:, ID) = 0;
    deIDs(nde + 1 : nde + 2) = elemChildren;
    nde = nde + 2;    
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

% Refresh face-element connections
msh.faceElems(1, remFIDs) = neIDs(msh.faceElems(1, remFIDs));
remFIDs_int = remFIDs(msh.faceType(remFIDs) == 0 | msh.faceType(remFIDs) == 1);
msh.faceElems(2, remFIDs_int) = neIDs(msh.faceElems(2, remFIDs_int));

% Finally remove the deleted elements and faces from the mesh
msh.elemCenter(:, deIDs)   = [];
msh.elemLength(:, deIDs)   = [];
msh.elemSize(deIDs)        = [];
msh.elemFaces(:, deIDs)    = [];
msh.elemLevel(deIDs)       = [];
msh.elemParent(deIDs)      = [];
msh.elemChildren(:, deIDs) = [];
msh.elemJac(:, deIDs)      = [];
msh.faceNormalx(:, dfIDs)  = [];
msh.faceType(dfIDs)        = [];
msh.faceElems(:, dfIDs)    = [];
msh.faceNums(:, dfIDs)     = [];
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
msh.LFaces  = 1 : msh.nFaces;
msh.nLFaces = msh.nFaces;

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






