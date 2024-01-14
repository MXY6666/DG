% Copyright (C) Xiangyi Meng.

function plotCuboidMesh(msh, pIDs, color, drawAxis)

if (nargin < 1)
    error('Not enough arguments')
end
if (msh.type ~= 301) && (msh.type ~= 302)
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

ct = msh.elemCenter(:, pIDs);
hx = msh.elemLength(1, pIDs);
hy = msh.elemLength(2, pIDs);
hz = msh.elemLength(3, pIDs);

X = [ct(1, :)-hx/2; ct(1, :)+hx/2; ct(1, :)+hx/2; ct(1, :)-hx/2; ct(1, :)-hx/2; ct(1, :)-hx/2; ct(1, :)+hx/2; ct(1, :)+hx/2; ct(1, :)-hx/2; ct(1, :)-hx/2];
Y = [ct(2, :)-hy/2; ct(2, :)-hy/2; ct(2, :)+hy/2; ct(2, :)+hy/2; ct(2, :)-hy/2; ct(2, :)-hy/2; ct(2, :)-hy/2; ct(2, :)+hy/2; ct(2, :)+hy/2; ct(2, :)-hy/2];
Z = [ct(3, :)-hz/2; ct(3, :)-hz/2; ct(3, :)-hz/2; ct(3, :)-hz/2; ct(3, :)-hz/2; ct(3, :)+hz/2; ct(3, :)+hz/2; ct(3, :)+hz/2; ct(3, :)+hz/2; ct(3, :)+hz/2];    
plot3(X, Y, Z, '--', 'Color', color, 'LineWidth', 1); hold on
X = [ct(1, :)+hx/2; ct(1, :)+hx/2];
Y = [ct(2, :)-hy/2; ct(2, :)-hy/2];
Z = [ct(3, :)-hz/2; ct(3, :)+hz/2];
plot3(X, Y, Z, '--', 'Color', color, 'LineWidth', 1); hold on
X = [ct(1, :)+hx/2; ct(1, :)+hx/2];
Y = [ct(2, :)+hy/2; ct(2, :)+hy/2];
Z = [ct(3, :)-hz/2; ct(3, :)+hz/2];
plot3(X, Y, Z, '--', 'Color', color, 'LineWidth', 1); hold on    
X = [ct(1, :)-hx/2; ct(1, :)-hx/2];
Y = [ct(2, :)+hy/2; ct(2, :)+hy/2];
Z = [ct(3, :)-hz/2; ct(3, :)+hz/2];
plot3(X, Y, Z, '--', 'Color', color, 'LineWidth', 1); hold off    

if ~drawAxis
    axis off
end



end





