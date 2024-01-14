% Copyright (C) Xiangyi Meng.

function plotLineMesh(msh, pIDs, nLabs)

if (nargin < 1)
    error('Not enough arguments')
end
if (msh.type ~= 101) && (msh.type ~= 102)
    error('Wrong mesh type')
end

% Default we plot all the leaf elements
if (nargin < 2) || isempty(pIDs)
    pIDs = msh.LElems;
end

if (nargin < 3) || isempty(nLabs)
    nLabs = 2;
end

ct = msh.elemCenter(:, pIDs);
h  = msh.elemLength(:, pIDs);
x  = unique([ct - h / 2, ct + h / 2]);

nx = numel(x);
if (nLabs > nx)
    nLabs = nx;
end
labs = cell(1, nx);
pos = round(linspace(1, nx, nLabs));
labs(pos) = num2cell(x(pos));

c = get(gcf,'Color');
ax = axes('YColor', c, 'Color', c, 'YTick', [], 'XTick', x);
ax.XLim = msh.dm;
ax.XTick = x;
ax.XTickLabel = labs;

end
