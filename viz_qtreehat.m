function [wp, xg, yg] = viz_qtreehat(Nx, Ny, srctype, potstr)
%
% function [wp, xg, yg] = viz_qtreehat(Nx, Ny, srctype, potstr)
%

if nargin == 0 
  Nx = 301;
  Ny = 301;
  srctype = 'impulse';
  potstr = 'log';
end

ep = 0.10;
epk = 0.70 * 2.0 / 300;
maxleaf = 8;

Xhr = 2.0;  % half range
xg = linspace(-Xhr, Xhr, Nx);
xg = xg(:);

Yhr = 2.0;
yg = linspace(-Yhr, Yhr, Ny);
yg = yg(:);

Xg = repmat(xg', [Ny, 1]);
Yg = repmat(yg, [1, Nx]);

xp = Xg(:);
yp = Yg(:);

if strcmp(srctype, 'impulse')
  x0 = 0.5;
  y0 = -0.5;
  vp = zeros(size(xp));
  [~, idx] = min((xp - x0).^2 + (yp - y0).^2);
  vp(idx) = 1.0;
elseif strcmp(srctype, 'randn')
  vp = randn(size(xp));
elseif strcmp(srctype, 'rand')
  vp = rand(size(xp));
elseif strcmp(srctype, 'rand-centered')
  vp = 2.0 * rand(size(xp)) - 1.0;
else
  error('unrecognized srctype');
end

wp = qtreehat(xp, yp, vp, ep, maxleaf, epk, potstr);

Wvalue = reshape(wp(:, 1), [Ny, Nx]);
Wgradx = reshape(wp(:, 2), [Ny, Nx]);
Wgrady = reshape(wp(:, 3), [Ny, Nx]);

if nargout == 0 
  figure;
  imagesc(xg, yg, Wgradx);
  colorbar;
  axis xy;
  xlabel('x', 'FontSize', 20);
  ylabel('y', 'FontSize', 20);
  title('gradx field', 'FontSize', 20);

  figure;
  imagesc(xg, yg, Wgrady);
  colorbar;
  axis xy;
  xlabel('x', 'FontSize', 20);
  ylabel('y', 'FontSize', 20);
  title('grady field', 'FontSize', 20);

  figure;
  imagesc(xg, yg, Wvalue);
  colorbar;
  axis xy;
  xlabel('x', 'FontSize', 20);
  ylabel('y', 'FontSize', 20);
  title('value field', 'FontSize', 20);
end

end
