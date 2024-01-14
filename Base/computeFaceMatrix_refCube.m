% leftFaceNum  : local face number of left element
% rightFaceNum : local face number of right element
% f2f          : face to face index
% flux         : numerical flux
% bs           : basis function set data (a struct) or just the degree of 
%                polynomial (an integer)
% mtol         : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function FM = computeFaceMatrix_refCube(leftFaceNum, rightFaceNum, f2f, flux, bs, mtol)

if (nargin < 4)
    error('Not enough arguments')
end

if (leftFaceNum <= 0)
    error('Wrong face number w.r.t. the left element')
end

if (rightFaceNum == 0)
    rightFaceNum = leftFaceNum;
end

% Default we choose basis functions of degree 1 
if (nargin < 5) || isempty(bs)
    bs = 1;
end

if (nargin < 6) || isempty(mtol)
    mtol = 1.0e-12;
end

if (leftFaceNum > 6) || (rightFaceNum > 6)
    isAdaptive = true;
else
    isAdaptive = false;
end

if isfloat(bs) && (bs >= 0)
    k     = floor(bs);
    quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);
    quad3 = GaussQuadratureRule_cube([k + 1, k + 1, k + 1], 301);
    bs    = setBasisFunctionSet_cube(quad2, quad3, 300 + k);
    if isAdaptive
        bs = addBasisFunctionSet_cube(bs, quad2, quad3);
    end
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'cube')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if ~(bs.facePointsType == 201 && 2 * bs.nfps(1) - 1 >= 2 * bs.deg && 2 * bs.nfps(2) - 1 >= 2 * bs.deg)   
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
        quad3 = GaussQuadratureRule_cube([bs.deg + 1, bs.deg + 1, bs.deg + 1], 301);
        bs    = setBasisFunctionSet_cube(quad2, quad3, bs.type);
        if isAdaptive
            bs = addBasisFunctionSet_cube(bs, quad2, quad3);
        end
    end    
else 
    error('Wrong argument bs')  
end

FM = cell(2, 2);
switch flux
    case 'l'
        FM{1, 1} =  bs.phitw_face{1, leftFaceNum}    * bs.phi_face{1, leftFaceNum};
        FM{2, 1} = -bs.phitw_face{f2f, rightFaceNum} * bs.phi_face{1, leftFaceNum};
        
        FM{1, 1} = mychop(FM{1, 1}, mtol);
        FM{2, 1} = mychop(FM{2, 1}, mtol);
    case 'r'
        FM{1, 2} =  bs.phitw_face{1, leftFaceNum}    * bs.phi_face{f2f, rightFaceNum};
        FM{2, 2} = -bs.phitw_face{f2f, rightFaceNum} * bs.phi_face{f2f, rightFaceNum};
        
        FM{1, 2} = mychop(FM{1, 2}, mtol);
        FM{2, 2} = mychop(FM{2, 2}, mtol);        
    case 'c'
        FM{1, 1} =  0.5 * bs.phitw_face{1, leftFaceNum}    * bs.phi_face{1, leftFaceNum};
        FM{1, 2} =  0.5 * bs.phitw_face{1, leftFaceNum}    * bs.phi_face{f2f, rightFaceNum};
        FM{2, 1} = -0.5 * bs.phitw_face{f2f, rightFaceNum} * bs.phi_face{1, leftFaceNum};
        FM{2, 2} = -0.5 * bs.phitw_face{f2f, rightFaceNum} * bs.phi_face{f2f, rightFaceNum}; 
        
        FM{1, 1} = mychop(FM{1, 1}, mtol);
        FM{1, 2} = mychop(FM{1, 2}, mtol);
        FM{2, 1} = mychop(FM{2, 1}, mtol);   
        FM{2, 2} = mychop(FM{2, 2}, mtol); 
    otherwise
        error('Wrong flux')
end


end





