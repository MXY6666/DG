% Copyright (C) Xiangyi Meng.

function plotRectMesh(msh, pIDs, color, ptype, showElemLabels, showFaceLabels, showFaceNormal, drawAxis)

if (nargin < 1)
    error('Not enough arguments')
end
if (msh.type ~= 201) && (msh.type ~= 202)
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

if (nargin < 7)  || isempty(showFaceNormal)
    showFaceNormal = false;
end

if (nargin < 8)  || isempty(drawAxis)
    drawAxis = true;
end

npe = length(pIDs);
ct = msh.elemCenter(:, pIDs);
hx = msh.elemLength(1, pIDs);
hy = msh.elemLength(2, pIDs);
X = [ct(1, :)-hx/2; ct(1, :)+hx/2; ct(1, :)+hx/2; ct(1, :)-hx/2; ct(1, :)-hx/2];
Y = [ct(2, :)-hy/2; ct(2, :)-hy/2; ct(2, :)+hy/2; ct(2, :)+hy/2; ct(2, :)-hy/2];
if (ptype == 1)
    plot(X, Y, color(1), 'LineWidth', 1);
else
    if (msh.type == 201)
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
    ect = [ct; zeros(1, npe)];
    elabs = 'e' + string(pIDs(:)); elabs = cellstr(elabs);
    et = matlab.graphics.primitive.world.Text('VertexData', single(ect), 'String', elabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    et.Font.Name = 'Helvetica';
    et.Font.Size = 10;
end

if showFaceLabels
    faceIDs = msh.elemFaces(:, pIDs);
    faceIDs = unique(faceIDs(:));
    if (msh.type == 202)
        faceIDs = faceIDs(msh.faceChildren(1, faceIDs) == 0);
    end
    leIDs = msh.faceElems(1, faceIDs);
    
    fct = [msh.elemCenter(:, leIDs); zeros(1, length(faceIDs))];
    IDs = find(msh.faceNums(1, faceIDs) == 1);
    fct(1, IDs) = fct(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 2;
    IDs = find(msh.faceNums(1, faceIDs) == 2);
    fct(1, IDs) = fct(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 2;  
    IDs = find(msh.faceNums(1, faceIDs) == 3);
    fct(2, IDs) = fct(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 2;    
    IDs = find(msh.faceNums(1, faceIDs) == 4);
    fct(2, IDs) = fct(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 2;     
    IDs = find(msh.faceNums(1, faceIDs) == 5);
    fct(1, IDs) = fct(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 2;
    fct(2, IDs) = fct(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 4;
    IDs = find(msh.faceNums(1, faceIDs) == 6);
    fct(1, IDs) = fct(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 2;
    fct(2, IDs) = fct(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 4;
    IDs = find(msh.faceNums(1, faceIDs) == 7);
    fct(1, IDs) = fct(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 2;
    fct(2, IDs) = fct(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 4; 
    IDs = find(msh.faceNums(1, faceIDs) == 8);
    fct(1, IDs) = fct(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 2;
    fct(2, IDs) = fct(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 4; 
    IDs = find(msh.faceNums(1, faceIDs) == 9);
    fct(1, IDs) = fct(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 4;
    fct(2, IDs) = fct(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 2;
    IDs = find(msh.faceNums(1, faceIDs) == 10);
    fct(1, IDs) = fct(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 4;
    fct(2, IDs) = fct(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 2;     
    IDs = find(msh.faceNums(1, faceIDs) == 11);
    fct(1, IDs) = fct(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 4;
    fct(2, IDs) = fct(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 2; 
    IDs = find(msh.faceNums(1, faceIDs) == 12);
    fct(1, IDs) = fct(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 4;
    fct(2, IDs) = fct(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 2;    
    
    flabs = 'f' + string(faceIDs); flabs = cellstr(flabs);
    ft = matlab.graphics.primitive.world.Text('VertexData', single(fct), 'String', flabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    ft.Font.Name = 'Helvetica';
    ft.Font.Size = 10;
end

if showFaceNormal
    faceIDs = msh.elemFaces(:, pIDs);
    faceIDs = unique(faceIDs(:));
    if (msh.type == 202)
        faceIDs = faceIDs(msh.faceChildren(1, faceIDs) == 0);
    end
    leIDs = msh.faces(1, faceIDs);
    
    xy = msh.elemCenter(:, leIDs);   
    IDs = find(msh.faceNums(1, faceIDs) == 1);
    xy(1, IDs) = xy(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 2;
    IDs = find(msh.faceNums(1, faceIDs) == 2);
    xy(1, IDs) = xy(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 2;   
    IDs = find(msh.faceNums(1, faceIDs) == 3);
    xy(2, IDs) = xy(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 2;    
    IDs = find(msh.faceNums(1, faceIDs) == 4);
    xy(2, IDs) = xy(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 2;    
    IDs = find(msh.faceNums(1, faceIDs) == 5);
    xy(1, IDs) = xy(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 2;
    xy(2, IDs) = xy(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 4;
    IDs = find(msh.faceNums(1, faceIDs) == 6);
    xy(1, IDs) = xy(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 2;
    xy(2, IDs) = xy(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 4;
    IDs = find(msh.faceNums(1, faceIDs) == 7);
    xy(1, IDs) = xy(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 2;
    xy(2, IDs) = xy(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 4; 
    IDs = find(msh.faceNums(1, faceIDs) == 8);
    xy(1, IDs) = xy(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 2;
    xy(2, IDs) = xy(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 4;  
    IDs = find(msh.faceNums(1, faceIDs) == 9);
    xy(1, IDs) = xy(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 4;
    xy(2, IDs) = xy(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 2;
    IDs = find(msh.faceNums(1, faceIDs) == 10);
    xy(1, IDs) = xy(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 4;
    xy(2, IDs) = xy(2, IDs) - msh.elemLength(2, leIDs(IDs)) / 2;     
    IDs = find(msh.faceNums(1, faceIDs) == 11);
    xy(1, IDs) = xy(1, IDs) - msh.elemLength(1, leIDs(IDs)) / 4;
    xy(2, IDs) = xy(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 2; 
    IDs = find(msh.faceNums(1, faceIDs) == 12);
    xy(1, IDs) = xy(1, IDs) + msh.elemLength(1, leIDs(IDs)) / 4;
    xy(2, IDs) = xy(2, IDs) + msh.elemLength(2, leIDs(IDs)) / 2;        
    
    nx = msh.faceNormalx(:, faceIDs);
    ny = msh.faceNormaly(:, faceIDs);

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
