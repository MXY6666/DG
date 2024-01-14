% leftFaceNum  : local face number of left element
% rightFaceNum : local face number of right element
% flux         : numerical flux
% bs           : basis function set data (a struct) or just the degree of 
%                polynomial (an integer)
% mtol         : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function FM = computeFaceMatrix_refLine(leftFaceNum, rightFaceNum, flux, bs, mtol)

if (nargin < 3)
    error('Not enough arguments')
end

if (leftFaceNum <= 0)
    error('Wrong face number w.r.t. the left element')
end

if (rightFaceNum == 0)
    rightFaceNum = leftFaceNum;
end

% Default we choose basis functions of degree 1 
if (nargin < 4) || isempty(bs)
    bs = 1;
end

if (nargin < 5) || isempty(mtol)
    mtol = 1.0e-12;
end

if isfloat(bs) && (bs >= 1)
    k    = floor(bs);
    quad = GaussQuadratureRule_line(k + 1, 101);  
    bs   = setBasisFunctionSet_line(quad, 100 + k);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'line')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
else
    error('Wrong argument bs')   
end

FM = cell(2, 2);
switch flux
    case 'l'
        FM{1, 1} =  bs.phitw_face{leftFaceNum}  * bs.phi_face{leftFaceNum};
        FM{2, 1} = -bs.phitw_face{rightFaceNum} * bs.phi_face{leftFaceNum};
        
        FM{1, 1} = mychop(FM{1, 1}, mtol);
        FM{2, 1} = mychop(FM{2, 1}, mtol);
    case 'r'
        FM{1, 2} =  bs.phitw_face{leftFaceNum}  * bs.phi_face{rightFaceNum};
        FM{2, 2} = -bs.phitw_face{rightFaceNum} * bs.phi_face{rightFaceNum};
        
        FM{1, 2} = mychop(FM{1, 2}, mtol);
        FM{2, 2} = mychop(FM{2, 2}, mtol);        
    case 'c'
        FM{1, 1} =  0.5 * bs.phitw_face{leftFaceNum}  * bs.phi_face{leftFaceNum};
        FM{1, 2} =  0.5 * bs.phitw_face{leftFaceNum}  * bs.phi_face{rightFaceNum};
        FM{2, 1} = -0.5 * bs.phitw_face{rightFaceNum} * bs.phi_face{leftFaceNum};
        FM{2, 2} = -0.5 * bs.phitw_face{rightFaceNum} * bs.phi_face{rightFaceNum}; 
        
        FM{1, 1} = mychop(FM{1, 1}, mtol);
        FM{1, 2} = mychop(FM{1, 2}, mtol);
        FM{2, 1} = mychop(FM{2, 1}, mtol);   
        FM{2, 2} = mychop(FM{2, 2}, mtol); 
    otherwise
        error('Wrong flux')
end


