function errs = line_qtreehat(potstr, N, maxleaf, NL, ep, epk)
%
% function errs = line_qtreehat(potstr, N, maxleaf, NL, ep, epk)
%
% Plot the result of line integrating the potential delta via its gradient
% compared to its direct evaluation, for a specific accuracy parameter "ep".
%
% N = number of source points, NL = number of line integration points
%
% EXAMPLE:
%   line_qtreehat('log', 100, 1, 1e3, 0.01, 0.001);
%

XY = randn(N, 2); % sources
Q = randn(N, 1);

V = randn(2, 1);
Vhat = V / norm(V);
s = linspace(-1.0, 1.0, NL);
s = s(:);
XYt = [s * Vhat(1), s * Vhat(2)];  % targets

wt_ref = qtreehat(XY(:, 1), XY(:, 2), Q, maxleaf, 0, epk, potstr, XYt(:, 1), XYt(:, 2));
pot_exact = wt_ref(:, 1);

wt = qtreehat(XY(:, 1), XY(:, 2), Q, maxleaf, ep, epk, potstr, XYt(:, 1), XYt(:, 2));
pot_direct = wt(:, 1);

fs = wt(:, 2) * Vhat(1) + wt(:, 3) * Vhat(2);
pot_delta = cumtrapz(s, fs);
pot_via_line = pot_direct(1) + pot_delta;

errs = [max(abs(pot_direct - pot_exact)), max(abs(pot_direct - pot_via_line))] / mean(abs(pot_exact));

figure;
plot(s, [pot_exact, pot_direct, pot_via_line], 'LineWidth', 3);
grid on;
hl = legend(sprintf('exact (%s)', potstr), sprintf('direct (ep=%.4f)', ep), 'line integral');
set(hl, 'FontSize', 15);
xlabel('line parameter s', 'FontSize', 20);
ylabel('potential value at (x(s), y(s))', 'FontSize', 20);
title(sprintf('N=%i sources, NL=%i line integral targets', N, NL), 'FontSize', 15);

end
