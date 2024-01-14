% Copyright (C) Xiangyi Meng.

function plotQuadMesh(msh, pIDs, color, ptype, showElemLabels, showFaceLabels, showNodeLabels, showFaceNormal, drawAxis)

if (nargin < 1)
    error('Not enough arguments')
end
if (msh.type ~= 211) && (msh.type ~= 212)
    error('Wrong mesh type')
end

% Default we plot all the leaf elements
if (nargin < 2) || isempty(pIDs)
    pIDs = msh.LElems;
end

if (nargin < 3)  || isempty(color)
    color = ['w', 'y', 'c', 'g', 'b', 'm', 'r', 'k'];
else
    fullColor = ['w', 'y', 'c', 'g', 'b', 'm', 'r', 'k'];
    ind1 = zeros(1, length(color));
    count = 0;
    for i = 1 : length(color)
        pos = strfind(fullColor, color(i));
        if ~isempty(pos)
            count = count + 1;
            ind1(count) = pos;
        end
    end
    ind1(count + 1 : end) = [];
    color = fullColor(ind1);
    
    ind2 = 1 : 8;
    ind2(ind1) = [];
    color(length(ind1) + 1 : 8) = fullColor(ind2);
end

if (nargin < 4)  || isempty(ptype)
    ptype = 1;
end
if (ptype ~= 1) && (ptype ~= 2)
    error('Wrong plot type')
end
if (ptype == 1) && (color(1) == 'w')
    color(1) = 'k';
end

if (nargin < 5)  || isempty(showElemLabels)
    showElemLabels = false;
end

if (nargin < 6)  || isempty(showFaceLabels)
    showFaceLabels = false;
end

if (nargin < 7)  || isempty(showNodeLabels)
    showNodeLabels = false;
end

if (nargin < 8)  || isempty(showFaceNormal)
    showFaceNormal = false;
end

if (nargin < 9)  || isempty(drawAxis)
    drawAxis = true;
end

npe = length(pIDs);
if (ptype == 1)
    Q = sparse(msh.elemNodes([1, 2, 3, 4], pIDs), msh.elemNodes([2, 3, 4, 1], pIDs), 1, msh.nNodes, msh.nNodes);
    [I, J]=find(Q | Q');
    K = find(I >= J);
    I = I(K);
    J = J(K);
    X = [msh.nodeCoor(1, I); msh.nodeCoor(1, J); NaN(1, length(I))];
    Y = [msh.nodeCoor(2, I); msh.nodeCoor(2, J); NaN(1, length(I))];
    X = X(:);
    Y = Y(:);

    plot(X, Y, color(1), 'LineWidth', 1);
else
    X = [msh.nodeCoor(1, msh.elemNodes(1, pIDs)); msh.nodeCoor(1, msh.elemNodes(2, pIDs)); msh.nodeCoor(1, msh.elemNodes(3, pIDs)); msh.nodeCoor(1, msh.elemNodes(4, pIDs))];
    Y = [msh.nodeCoor(2, msh.elemNodes(1, pIDs)); msh.nodeCoor(2, msh.elemNodes(2, pIDs)); msh.nodeCoor(2, msh.elemNodes(3, pIDs)); msh.nodeCoor(2, msh.elemNodes(4, pIDs))];
    if (msh.type == 211)
        patch(X, Y, color(1), 'LineWidth', 1);
    else
        level = msh.elemLevel(pIDs);
        rgb = ones(npe, 3);
        for i = 0 : msh.maxLevel
            pos = find(level == i);
            if (i < 8)
                rgb(pos, :) = repmat(c2rgb(color(i + 1)), [length(pos), 1]);
            end
        end
        rgb = reshape(rgb, [npe, 1, 3]);
        patch(X, Y, rgb, 'LineWidth', 1);                
    end
end


if showElemLabels
    ect = [0.25 * (msh.nodeCoor(:, msh.elemNodes(1, pIDs)) + msh.nodeCoor(:, msh.elemNodes(2, pIDs)) + msh.nodeCoor(:, msh.elemNodes(3, pIDs)) + msh.nodeCoor(:, msh.elemNodes(4, pIDs))); zeros(1, npe)];
    elabs = 'e' + string(pIDs(:)); elabs = cellstr(elabs);
    et = matlab.graphics.primitive.world.Text('VertexData', single(ect), 'String', elabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    et.Font.Name = 'Helvetica';
    et.Font.Size = 10;
end

if showFaceLabels
    faceIDs = msh.elemFaces(:, pIDs);
    faceIDs = unique(faceIDs(:));
    if (msh.type == 212)
        faceIDs = faceIDs(msh.faceChildren(1, faceIDs) == 0);
    end
    fct = [0.5 * (msh.nodeCoor(:, msh.faceNodes(1, faceIDs)) + msh.nodeCoor(:, msh.faceNodes(2, faceIDs))); zeros(1, length(faceIDs))];   
    flabs = 'f' + string(faceIDs); flabs = cellstr(flabs);
    ft = matlab.graphics.primitive.world.Text('VertexData', single(fct), 'String', flabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    ft.Font.Name = 'Helvetica';
    ft.Font.Size = 10;
end

if showNodeLabels   
    nodeIDs = msh.elemNodes(:, pIDs);
    nodeIDs = unique(nodeIDs(:));
    nc = [msh.nodeCoor(:, nodeIDs); zeros(1, length(nodeIDs))];
    nlabs = 'n' + string(nodeIDs); nlabs = cellstr(nlabs);
    nt = matlab.graphics.primitive.world.Text('VertexData', single(nc), 'String', nlabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    nt.Font.Name = 'Helvetica';
    nt.Font.Size = 10;    
end

if showFaceNormal
    faceIDs = msh.elemFaces(:, pIDs);
    faceIDs = unique(faceIDs(:));
    if (msh.type == 212)
        faceIDs = faceIDs(msh.faceChildren(1, faceIDs) == 0);
    end
    nodes1 = msh.nodeCoor(:, msh.faceNodes(1, faceIDs));
    nodes2 = msh.nodeCoor(:, msh.faceNodes(2, faceIDs));
    xy = 0.5 * (nodes1 + nodes2);
    h = sqrt(sum((nodes2 - nodes1).^2));
    nx = (nodes2(2, :) - nodes1(2, :)) ./ h;
    ny = (nodes1(1, :) - nodes2(1, :)) ./ h;

    hold on 
    quiver(xy(1, :), xy(2, :), nx, ny, 'r-');
    hold off
end

if ~drawAxis
    axis off
end

function rgb = c2rgb(c)
    switch c
        case 'r'
            rgb = [1 0 0];
        case 'g'
            rgb = [0 1 0];
        case 'b'
            rgb = [0 0 1];
        case 'c'
            rgb = [0 1 1];
        case 'm'
            rgb = [1 0 1];
        case 'y'
            rgb = [1 1 0];
        case 'k'
            rgb = [0 0 0];
        case 'w'
            rgb = [1 1 1];
        otherwise
            error('Wrong short name for color')
    end
end

end
