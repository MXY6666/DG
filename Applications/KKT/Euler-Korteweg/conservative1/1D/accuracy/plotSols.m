function plotSols(msh, U, t, tc, quad1, bs)

U = reshape(U, [bs.nb, 2 * msh.nLElems]);
% center of every element
ct  = msh.elemCenter(:, msh.LElems);

% exact solutions
rho = tc.rho(ct, t);
u   = tc.u(ct, t);

% mean values of numerical solutions
rhoh = 0.5 * (quad1.weights' * (bs.phi{1} * U(:, 1 : msh.nLElems)));
uh   = 0.5 * (quad1.weights' * (bs.phi{1} * U(:, msh.nLElems + 1 : end)));

% Plot
subplot(1, 2, 1)
plot(ct, rho, '-b', 'LineWidth', 1); hold on
plot(ct, rhoh, '*r')
xlabel('X', 'FontSize', 16);
ylabel('Density', 'FontSize', 16);
legend('\rho', '\rho_h', 'FontSize', 16);
hold off

subplot(1, 2, 2)
plot(ct, u, '-b', 'LineWidth', 1); hold on
plot(ct, uh, '*r')
xlabel('X', 'FontSize', 16)
ylabel('Velocity', 'FontSize', 16);
legend('u', 'u_h', 'FontSize', 16);
hold off

end