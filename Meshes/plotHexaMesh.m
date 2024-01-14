% Copyright (C) Xiangyi Meng.

function plotHexaMesh(msh, pIDs, color, drawAxis)

if (nargin < 1)
    error('Not enough arguments')
end
if (msh.type ~= 311) && (msh.type ~= 312)
    error('Wrong mesh type')
end

% Default we plot all the leaf elements
if (nargin < 2) || isempty(pIDs)
    pIDs = msh.LElems;
end

if (nargin < 3)  || isempty(color)
    color = 'k';
end

if (nargin < 4)  || isempty(drawAxis)
    drawAxis = true;
end

v1 = [1, 2, 3, 4, 1, 5, 6, 7, 8, 2, 3, 4];
v2 = [2, 3, 4, 1, 5, 6, 7, 8, 5, 6, 7, 8];
Q = sparse(msh.elemNodes(v1, pIDs), msh.elemNodes(v2, pIDs), 1, msh.nNodes, msh.nNodes);
[I, J] = find(Q | Q');
K = find(I >= J);
I = I(K);
J = J(K);
X = [msh.nodeCoor(1, I); msh.nodeCoor(1, J); NaN(1, length(I))];
Y = [msh.nodeCoor(2, I); msh.nodeCoor(2, J); NaN(1, length(I))];
Z = [msh.nodeCoor(3, I); msh.nodeCoor(3, J); NaN(1, length(I))];
X = X(:);
Y = Y(:);
Z = Z(:);
plot3(X, Y, Z, '--', 'Color', color, 'LineWidth', 1);

if ~drawAxis
    axis off
end


end