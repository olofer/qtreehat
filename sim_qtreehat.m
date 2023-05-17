function rep = sim_qtreehat(N, Nt, dt, ep, potstr)
%
% function rep = sim_qtreehat(N, Nt, dt, ep, potstr)
%
% 2D log-potential gravity simulation with N particles.
% 2D 1/r potential is also supported.
%
% Conservative simulation: E = K + V = constant.
%   K = 0.5 * sum_i m_i * v_i^2
%   V = sum_i m_i * phi_i
%
% Here the gravitational field phi (at ith particle) is calculated using qtreehat(.):
%   phi(r_i) = sum_j m_j * log |r_i - r_j|
%
% Or if potstr = 'inv':
%   phi(r_i) = -sum_j m_j / |r_i - r_j|
%
% The units are implicitly such that the factor in front of the field is 1.
%
% Equation of motion:
%   dot(r_i) = v_i
%   dot(v_i) = -grad phi(r_i)
%

%{
  
EXAMPLE:
  close all; sim_qtreehat(256, 10000, 1e-3, 0.15); rep=ans;

It seems that the treecode approximation degrades momentum conservation faster than energy conservation.
(can be an unreliable observation).

NOTE: there is some issue with accuracy of energy conservation when using 'inv' potential, not sure why..
  octave:3> sim_qtreehat(256, 10000, 1e-4, 0, 'inv');
  octave:4> sim_qtreehat(256, 10000, 1e-4, 0, 'log');

%}

if nargin < 5
  potstr = 'log';
end

rep = struct;
rep.creator = mfilename();
rep.energy = NaN(Nt, 2);
rep.linmom = NaN(Nt, 2);
rep.angmom = NaN(Nt, 1);
rep.time = NaN(Nt, 1);
rep.params = struct;
rep.params.ep = ep;
rep.params.maxleaf = 8;
rep.params.M = 7.35e22;
rep.params.L = 3.736e8;
rep.params.epk = 1.73e6 / rep.params.L;    % radius of moon / earth-moon distance
rep.params.G = 6.674e-11;  % https://en.wikipedia.org/wiki/Gravitational_constant
rep.params.potstr = potstr;

% dim of C = pi * params.G * params.M * parms.T^2 / params.L ^ 2;
% is [1/L] (per meter) which I think is what it should be
% so there is a time-scale T where this coefficient is exactly 1.
% this timescale is the time unit in the below simulation.

rep.params.T = sqrt(rep.params.L^2 / (pi * rep.params.G * rep.params.M));

Mmax = 1.0;
Mmin = 0.1 * Mmax;
Xrange = 1.0; % distance earth-moon
Yrange = Xrange;

M = Mmin + (Mmax - Mmin) * rand(N, 1); 
X = (rand(N, 1) - 0.5) * Xrange;
Y = (rand(N, 1) - 0.5) * Yrange;
[Ax, Ay, Phi] = calc_acceleration(X, Y, M, rep.params);
VX = zeros(N, 1);
VY = zeros(N, 1);

t = 0;
for k = 1:Nt
  [K, V] = calc_energy(M, VX, VY, Phi);
  [PX, PY] = calc_momentum(VX, VY, M);
  LZ = calc_ang_momentum(X, VX, Y, VY, M);
  rep.energy(k, :) = [K, V];
  rep.linmom(k, :) = [PX, PY];
  rep.angmom(k) = LZ;
  rep.time(k) = t;

  % https://en.wikipedia.org/wiki/Leapfrog_integration

  VXh = VX + dt * Ax / 2;
  VYh = VY + dt * Ay / 2;
  X = X + dt * VXh;
  Y = Y + dt * VYh;
  [Ax, Ay, Phi] = calc_acceleration(X, Y, M, rep.params);
  VX = VXh + dt * Ax / 2;
  VY = VYh + dt * Ay / 2;
  t = t + dt;

  if rem(k, 100) == 0 && nargout == 0 && N < 100 && ep == 0
    figure(1); clf;
    plot(X, Y, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10);
    axis equal;
    title(['time = ', num2str(t)], 'FontSize', 20);
    pause(1);
  end

end

if nargout == 0 

  mtot = sum(M);
  tsec = rep.params.T * rep.time;

  figure;
  plot(tsec, [rep.linmom, rep.angmom] / mtot, 'LineWidth', 3);
  hl = legend('x', 'y', 'angular');
  set(hl, 'FontSize', 20);
  xlabel('time [seconds]', 'FontSize', 20);
  ylabel('total momentum / total mass', 'FontSize', 20);
  title(sprintf('N=%i, time-unit=%.4f, Nt=%i, potstr=%s', N, rep.params.T, Nt, rep.params.potstr), 'FontSize', 20);

  tot = sum(rep.energy, 2);
  avgTot = mean(tot);

  figure;
  plot(tsec, [rep.energy, tot], 'LineWidth', 3);
  hl = legend('kinetic', 'potential', 'sum');
  set(hl, 'FontSize', 20);
  xlabel('time [seconds]', 'FontSize', 20);
  ylabel('normalized total energy', 'FontSize', 20);
  title(sprintf('N=%i, time-unit=%.4f, Nt=%i, potstr=%s', N, rep.params.T, Nt, rep.params.potstr), 'FontSize', 20);
end

end

function [Ax, Ay, Phi] = calc_acceleration(X, Y, M, params)
  W = qtreehat(X, Y, M, params.ep, params.maxleaf, params.epk, params.potstr);
  Ax = -1.0 * W(:, 2);
  Ay = -1.0 * W(:, 3);
  Phi = W(:, 1);
end

function [K, V] = calc_energy(M, VX, VY, Phi)
  K = sum(M .* (VX.^2 + VY.^2)) / 2;
  V = sum(M .* Phi) / 2;
end

function [PX, PY] = calc_momentum(VX, VY, M)
  PX = sum(M .* VX);
  PY = sum(M .* VY);
end

function LZ = calc_ang_momentum(X, VX, Y, VY, M)
  LZ = sum(M .* (X .* VY - Y .* VX));
end
