function mshPlt(x, y)

% number of points in xi and eta direction
[m, n] = size(x);

% Plot 
X(m * n + n + 1 : 2 * m * n + m + n) = [x, NaN(m, 1)]';
Y(m * n + n + 1 : 2 * m * n + m + n) = [y, NaN(m, 1)]';
X(1 : m * n + n) = [x; NaN(1, n)];
Y(1 : m * n + n) = [y; NaN(1, n)];

plot(X, Y, 'k', 'LineWidth', 1);

end