% the 9th row of NEI stores the number of neighbors
% Copyright (C) Xiangyi Meng.

function NEI = getElementNeighbor_rect(msh, elementIDs)

if (msh.type ~= 201) && (msh.type ~= 202) 
    error('wrong mesh type')
end

if (nargin < 2)  || isempty(elementIDs)
    elementIDs = msh.LElems;
end

ne = length(elementIDs);
NEI = zeros(9, ne);
neiIDs = zeros(9, 1);
for i = 1 : ne 
    ID = elementIDs(i);
    elemFaces    = msh.elemFaces(:, ID);
    faceTypes    = msh.faceType(:, elemFaces);
    faceElems    = msh.faceElems(:, elemFaces);
    faceChildren = msh.faceChildren(:, elemFaces);
        
    neiIDs(:) = 0;
    nnei = 0;
    % the left face
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        if (faceChildren(1, 1) ~= 0)
            neiIDs(nnei + 1 : nnei + 2) = msh.faceElems(1, faceChildren(:, 1));
            nnei = nnei + 2;
        else
            neiIDs(nnei + 1) = faceElems(1, 1);
            nnei = nnei + 1;            
        end        
    end 
    % the right face
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        if (faceChildren(1, 2) ~= 0)
            neiIDs(nnei + 1 : nnei + 2) = msh.faceElems(2, faceChildren(:, 2));
            nnei = nnei + 2;
        else
            neiIDs(nnei + 1) = faceElems(2, 2);
            nnei = nnei + 1;            
        end        
    end   
    % the bottom face
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        if (faceChildren(1, 3) ~= 0)
            neiIDs(nnei + 1 : nnei + 2) = msh.faceElems(1, faceChildren(:, 3));
            nnei = nnei + 2;
        else
            neiIDs(nnei + 1) = faceElems(1, 3);
            nnei = nnei + 1;
        end
    end    
    % the top face
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        if (faceChildren(1, 4) ~= 0)
            neiIDs(nnei + 1 : nnei + 2) = msh.faceElems(2, faceChildren(:, 4));
            nnei = nnei + 2;
        else
            neiIDs(nnei + 1) = faceElems(2, 4);
            nnei = nnei + 1;              
        end        
    end
    neiIDs(9) = nnei;
    NEI(:, i) = neiIDs;
end


end