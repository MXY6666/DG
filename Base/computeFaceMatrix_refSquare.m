% leftFaceNum  : local face number of left element
% rightFaceNum : local face number of right element
% f2f          : face to face index
% flux         : numerical flux
% bs           : basis function set data (a struct) or just the degree of 
%                polynomial (an integer)
% mtol         : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function FM = computeFaceMatrix_refSquare(leftFaceNum, rightFaceNum, f2f, flux, bs, mtol)

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

if (leftFaceNum > 4) || (rightFaceNum > 4)
    isAdaptive = true;
else
    isAdaptive = false;
end

if isfloat(bs) && (bs >= 0)
    k     = floor(bs);
    quad1 = GaussQuadratureRule_line(k + 1, 101);
    quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);
    bs    = setBasisFunctionSet_square(quad1, quad2, 200 + k);
    if isAdaptive
        bs = addBasisFunctionSet_square(bs, quad1, quad2);
    end
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'square')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if ~(bs.facePointsType == 101 && 2 * bs.nfp - 1 >= 2 * bs.deg) && ...
       ~(bs.facePointsType == 102 && 2 * bs.nfp - 3 >= 2 * bs.deg)     

        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type);
        if isAdaptive
            bs = addBasisFunctionSet_square(bs, quad1, quad2);
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



