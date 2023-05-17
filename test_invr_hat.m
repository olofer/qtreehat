function err_table = test_invr_hat(npt, farness)
%
% function err_table = test_invr_hat(npt, farness)
%
% Verify quadrupole expansion of a set of sources with (-)1/r potential.
% Plots are only generated if no output is requested explicitly.
%

if nargin < 2
  farness = 10.0;
end

assert(farness > 1.0);

x = randn(npt, 1);
y = randn(npt, 1);
q = rand(npt, 1); 

mx = mean(x);
my = mean(y);

rx = x - mx;
ry = y - my;

[qsum, qx, qy, qxx, qxy, qyy] = calc_summaries(q, rx, ry);

rmax = max(sqrt(sum([rx.^2, ry.^2], 2)));
Reval = farness * rmax;

th = linspace(-pi, pi, 200);
th = th(:);
Rx = cos(th) * Reval;
Ry = sin(th) * Reval;

[Vtrue, gxtrue, gytrue] = eval_potential(q, rx, ry, Rx, Ry);
[Vquad, gxquad, gyquad] = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, Rx, Ry);

err_table = NaN(3, 3);
for ordr = 0:2
  err_table(1 + ordr, 1) = relative_approx_error(Vtrue, Vquad, ordr);
  err_table(1 + ordr, 2) = relative_approx_error(gxtrue, gxquad, ordr);
  err_table(1 + ordr, 3) = relative_approx_error(gytrue, gyquad, ordr);
end

for c = 1:size(err_table, 2)
  if ~issorted(flipud(err_table(:, c)))
    warning('error table has acolumn which is not monotonic');
  end
end

fdep = 1e-4;
[~, gradX, gradY, fdX, fdY] = fdiff_eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, Rx, Ry, fdep);

fderrX = max(abs(fdX - gradX));
fderrY = max(abs(fdY - gradY));

if max(fderrX) > 1e-10 || max(fderrY) > 1e-10 
  warning('quadrupole expansion finit difference error detected');
end

if nargout > 0 
  return;
end

figure;
plot(x, y, 'k.');
hold on;
plot(mx + Rx, my + Ry, 'b-');
plot(mx + rmax*cos(th), my + rmax*sin(th), 'c-.');
xlabel('x');
ylabel('y');
grid on;
axis equal;
legend('source points', 'evaluation', 'rmax');

figure;
plot(th, [gxtrue, gytrue, ...
          gxquad(:, 1), gyquad(:, 1), ...
          sum(gxquad(:, 1:2), 2), sum(gyquad(:, 1:2), 2), ...
          sum(gxquad, 2), sum(gyquad, 2)], ...
     'LineWidth', 3);
xlabel('theta');
ylabel('Potential');
legend('gX(true)', 'gY(true)', ...
       'gX(0th)', 'gY(0th)', ...
       'gX(1st))', 'gY(1st)', ...
       'gX(2nd)', 'gY(2nd)');
title(sprintf('gradients @ rmax/Reval = %.4f', rmax / Reval));

figure;
plot(th, [Vtrue, Vquad(:, 1), sum(Vquad(:, 1:2), 2), sum(Vquad, 2)], 'LineWidth', 3);
xlabel('theta');
ylabel('Potential');
legend('Vtrue', 'Single (0th)', 'Dipole (1st)', 'Quadrupole (2nd)');
title(sprintf('potential @ rmax/Reval = %.4f; dip. err. = %.4e, quad. rerr. = %.4e', ...
              rmax / Reval, err_table(1, 2), err_table(1, 3)));

end

function [V, dVx, dVy] = eval_potential(q, rx, ry, Rx, Ry)
  S = numel(q);
  assert(numel(rx) == S && numel(ry) == S);
  T = numel(Rx);
  assert(numel(Ry) == T);
  V = zeros(T, 1);
  dVx = zeros(T, 1);
  dVy = zeros(T, 1);
  for ii = 1:S
    Dsisq = sum([Rx - rx(ii), Ry - ry(ii)].^2, 2);
    Dsi = sqrt(Dsisq);
    V = V - q(ii) ./ Dsi;
    dVx = dVx + q(ii) * (Rx - rx(ii)) ./ (Dsisq .* Dsi);
    dVy = dVy + q(ii) * (Ry - ry(ii)) ./ (Dsisq .* Dsi);
  end
end

function [V, dVx, dVy] = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, Rx, Ry)
  T = numel(Rx);
  V = zeros(T, 3);
  dVx = zeros(T, 3);
  dVy = zeros(T, 3);
  for jj = 1:T
    Xj = Rx(jj);
    Yj = Ry(jj);
    Rj = sqrt(Xj^2 + Yj^2);
    V0 = -1 * qsum / Rj;
    V(jj, 1) = V0;

    Xhatj = Xj / Rj;
    Yhatj = Yj / Rj;
    V1 = -1 * (qx * Xhatj + qy * Yhatj) / Rj^2;
    V(jj, 2) = V1;

    Cxx = .5 * Yhatj^2 - Xhatj^2;
    Cyy = .5 * Xhatj^2 - Yhatj^2;
    Cxy = -3 * Xhatj * Yhatj;
    V2 = (Cxx * qxx + Cyy * qyy + Cxy * qxy) / Rj^3;
    V(jj, 3) = V2;

    dVx(jj, 1) = qsum * Xhatj / Rj^2;
    dVy(jj, 1) = qsum * Yhatj / Rj^2;

    CX_x = 2 * Xhatj^2 - Yhatj^2;
    CX_y = 3 * Xhatj * Yhatj;

    CY_x = 3 * Xhatj * Yhatj;
    CY_y = 2 * Yhatj^2 - Xhatj^2;

    dVx(jj, 2) = (qx * CX_x + qy * CX_y) / Rj^3;
    dVy(jj, 2) = (qx * CY_x + qy * CY_y) / Rj^3;

    CX_xx = 3 * Xhatj^3 - (9/2) * Xhatj * Yhatj^2;
    CY_xx = 6 * Xhatj^2 * Yhatj - (3/2) * Yhatj^3;

    CX_xy = -3 * Yhatj * (Yhatj^2 - 4 * Xhatj^2);
    CY_xy = -3 * Xhatj * (Xhatj^2 - 4 * Yhatj^2);

    CX_yy = 6 * Xhatj * Yhatj^2 - (3/2) * Xhatj^3;
    CY_yy = 3 * Yhatj^3 - (9/2) * Xhatj^2 * Yhatj;

    dVx(jj, 3) = (qxx * CX_xx + qxy * CX_xy + qyy * CX_yy) / Rj^4;
    dVy(jj, 3) = (qxx * CY_xx + qxy * CY_xy + qyy * CY_yy) / Rj^4;
  end
end

function [qsum, qx, qy, qxx, qxy, qyy] = calc_summaries(q, rx, ry)
  qsum = sum(q);
  qx = sum(q.*rx);
  qy = sum(q.*ry);
  qxx = sum(q.*rx.*rx);
  qxy = sum(q.*rx.*ry);
  qyy = sum(q.*ry.*ry);
end

% ordr should be 0, 1, or 2.
function err = relative_approx_error(V0, Vq, ordr)
  assert(ordr == 0 || ordr == 1 || ordr == 2);
  assert(size(V0, 2) == 1);
  assert(size(Vq, 2) == 3);
  assert(size(V0, 1) == size(Vq, 1));
  err = max(abs(V0 - sum(Vq(:, 1:(ordr + 1)), 2))) / max(abs(V0));
end

function [Value, gradX, gradY, fdX, fdY] = fdiff_eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, Rx, Ry, fdep) 
  [Value, gradX, gradY] = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, Rx, Ry);
  Npts = numel(Rx);
  fdX = NaN(Npts, 3);
  fdY = NaN(Npts, 3);
  for jj = 1:Npts
    VjxA = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, [Rx(jj) + fdep], [Ry(jj)]);
    VjxB = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, [Rx(jj) - fdep], [Ry(jj)]);
    fdX(jj, :) = (VjxA - VjxB) / (2 * fdep);

    VjyA = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, [Rx(jj)], [Ry(jj) + fdep]);
    VjyB = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, [Rx(jj)], [Ry(jj) - fdep]);
    fdY(jj, :) = (VjyA - VjyB) / (2 * fdep);
  end
end
