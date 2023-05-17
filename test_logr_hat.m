function err_table = test_logr_hat(npt, farness)
%
% function err_table = test_logr_hat(npt, farness)
%
% Numerical sanity check of expansion formulas for log(|r|) potential.
% Plots are only generated if no output is requested explicitly.
%

if nargin < 2
  farness = 10.0;
end

assert(farness > 1.0);

x = randn(npt, 1);
y = randn(npt, 1);
q = rand(npt, 1); % 2.0 * rand(npt, 1) - 1.0;   % + bias argument ?

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

if nargout == 1 
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
    Dsisq = sum([rx(ii) - Rx, ry(ii) - Ry].^2, 2);
    Dsi = sqrt(Dsisq);
    V = V + q(ii) * log(Dsi);
    dVx = dVx - q(ii) * (rx(ii) - Rx) ./ Dsisq;
    dVy = dVy - q(ii) * (ry(ii) - Ry) ./ Dsisq;
  end
end

function [V, dVx, dVy] = eval_quadrupole(qsum, qx, qy, qxx, qxy, qyy, Rx, Ry)
  T = numel(Rx);
  V = zeros(T, 3);  % add up columns to get quadrupole approx
  dVx = zeros(T, 3);
  dVy = zeros(T, 3);
  for jj = 1:T
    Xj = Rx(jj);
    Yj = Ry(jj);
    Rj = sqrt(Xj^2 + Yj^2);
    V0 = qsum * log(Rj);
    V(jj, 1) = V0;

    Xhatj = Xj / Rj;
    Yhatj = Yj / Rj;
    V1 = -1 * (qx * Xhatj + qy * Yhatj) / Rj;
    V(jj, 2) = V1;

    Cxx = -Xhatj^2 + Yhatj^2;
    Cyy = Xhatj^2 - Yhatj^2;
    Cxy = -2 * Xhatj * Yhatj;
    V2 = 0.5 * (Cxx * qxx + Cyy * qyy + Cxy * qxy) / Rj^2;
    V(jj, 3) = V2;

    dVx(jj, 1) = qsum * Xhatj / Rj;
    dVy(jj, 1) = qsum * Yhatj / Rj;

    dVx(jj, 2) = -1 * (qx * Cxx + qy * Cxy) / Rj^2;
    dVy(jj, 2) = -1 * (qx * Cxy + qy * Cyy) / Rj^2;

    dVx(jj, 3) = 0.5 * (2 * Xhatj * (Xhatj^2 - 3 * Yhatj^2) * (qxx - qyy) - 2 * Yhatj * (Yhatj^2 - 3 * Xhatj^2) * qxy) / Rj^3;
    dVy(jj, 3) = 0.5 * (2 * Yhatj * (Yhatj^2 - 3 * Xhatj^2) * (qyy - qxx) - 2 * Xhatj * (Xhatj^2 - 3 * Yhatj^2) * qxy) / Rj^3;
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
