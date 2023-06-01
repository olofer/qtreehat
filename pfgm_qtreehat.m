%
% Demonstration script illustrating crude (PFGM) Poisson Flow Generative Modeling for 1D data
% using 2D poisson equation solutions evaluated using the tree code MEX interface.
%

FSZ = 20;
NCT = 0; %15;

N = 2e3;

xp = zeros(N, 1);
yp = randn(N, 1);
%yp = -.5 + rand(N, 1);
vp = ones(N, 1) / N;

pot = 'log';
ep  = 0.10;
epk = 1e-4;
maxleaf = 32;

fprintf(1, '<ydata>   = %f\n', mean(yp));
fprintf(1, 'sd(ydata) = %f\n', std(yp));

% Draw M samples on a far-field circle of radius R
R = 10.0;
M = 500;
K = 5000;
ds = 2e-3;
theta = rand(M, 1) * 2 * pi;
XY = R * [cos(theta), sin(theta)];
P = NaN(M, 2, K);
fprintf(1, 'sampling %i flow paths (each for %i steps)..\n', M, K);
for k = 1:K
  P(:, :, k) = XY;
  wk = qtreehat(xp, yp, vp, maxleaf, ep, epk, pot, XY(:, 1), XY(:, 2));
  gx = wk(:, 2);
  gy = wk(:, 3);
  glen = sqrt(gx.^2 + gy.^2);
  ghat = [gx ./ glen, gy ./ glen];
  XY = XY - ghat * ds;
end

% Sample y-coordinate of final set of points as a 1D distribution
yfinal = squeeze(P(:, 2, K));
fprintf(1, '<yfinal>   = %f\n', mean(yfinal));
fprintf(1, 'sd(yfinal) = %f\n', std(yfinal));

ymin = min(yp);
ymax = max(yp);
yvec = linspace(ymin, ymax, 100).';
cdf_source = empirical_cdf(yvec, yp);
cdf_sample = empirical_cdf(yvec, yfinal);

figure;
plot(yvec, [cdf_source, cdf_sample], 'LineWidth', 3);
grid on;
hl = legend(sprintf('source (sz=%i)', N), ...
            sprintf('flow sample (sz=%i)', M));
set(hl, 'FontSize', FSZ);
xlabel('y', 'FontSize', FSZ);
ylabel('CDF', 'FontSize', FSZ);
title('Source CDF vs Flow sample CDF', 'FontSize', FSZ);

col = viridis(M);
figure;
hold on;
for m = 1:M 
  xym = squeeze(P(m, :, :));
  plot(xym(1, :), xym(2, :), 'Color', col(m, :), 'LineWidth', 2);
end
axis equal;
grid on;
xlabel('x', 'FontSize', FSZ);
ylabel('y', 'FontSize', FSZ);
title('Poisson Flow sample paths', 'FontSize', FSZ);

if NCT <= 0
  return;
end

nxt = 600;
nyt = 900;
xtvec = linspace(-4, 4, nxt).';
ytvec = linspace(-6, 6, nyt).';
Xt = repmat(xtvec.', [nyt, 1]);
Yt = repmat(ytvec, [1, nxt]);

qtreehat(xp, yp, vp, maxleaf, ep, epk, pot, Xt(:), Yt(:));
wt = ans;

phi = reshape(wt(:, 1), [nyt, nxt]);
phix = reshape(wt(:, 2), [nyt, nxt]);
phiy = reshape(wt(:, 3), [nyt, nxt]);

figure;
contourf(xtvec, ytvec, phi, NCT);
axis equal;
colorbar;
xlabel('x', 'FontSize', FSZ);
ylabel('y', 'FontSize', FSZ);
title('Contours of Poisson potential', 'FontSize', FSZ);
