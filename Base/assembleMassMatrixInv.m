% Assemble inverse of mass matrix
% Copyright (C) Xiangyi Meng.

function IM = assembleMassMatrixInv(msh, mtol, bs, IME)

t = round(mod(msh.type, 100) / 10);
switch t    
    case 0
        if (nargin < 2)
            bs = 1;
        end
        if (nargin < 3)
            IME = [];
        end
        h = round(msh.type / 100);
        switch h
            case 1
                IM = assembleElementMatrix_line(msh, 0, 0, bs, true, IME);
            case 2
                IM = assembleElementMatrix_rect(msh, [0; 0], [0; 0], bs, true, IME);
            case 3
                IM = assembleElementMatrix_cuboid(msh, [0; 0; 0], [0; 0; 0], bs, true, IME);
        end
    case 1
        nb = size(msh.massMatInv{1}, 1);
        re = repmat((1 : nb)', [nb, 1]);
        ce = repelem((1 : nb)', nb);
        
        r = re + (0 : msh.nLElems - 1) * nb;
        c = ce + (0 : msh.nLElems - 1) * nb;
        v = mychop([msh.massMatInv{msh.LElems}], mtol);
        IM = sparse(r, c, v, msh.nLElems * nb, msh.nLElems * nb);
end

end