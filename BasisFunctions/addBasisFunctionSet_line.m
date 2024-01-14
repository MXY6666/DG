% Evaluate basis functions at Gauss points of two subelements of reference 
% line
%            o------o------o
%                1      2
% Copyright (C) Xiangyi Meng.

function bs = addBasisFunctionSet_line(bs, quad1)

bs.phi_sube = {basisFunctionSet_line(-0.5 + 0.5 * quad1.points, bs.type, 0), ...
               basisFunctionSet_line( 0.5 + 0.5 * quad1.points, bs.type, 0)};
           
           
bs.phitw_sube = cell(1, 2);
for j = 1 : 2
    bs.phitw_sube{j} = (quad1.weights .* bs.phi_sube{j})';
end

end