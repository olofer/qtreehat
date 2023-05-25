function [epvec, eperr] = test_qtreehat(potstr, N, maxleaf, epk)
%
% function [epvec, eperr] = test_qtreehat(potstr, N, maxleaf, epk)
%
% Basic test script for the qtreehat.c MEX program.
% should show that results converge to reference for small ep,
% but is much faster for larger ep while still accurate to several digits
%
% EXAMPLES:
%   test_qtreehat('log', 8e3, 4);
%   test_qtreehat('inv', 8e3, 4);
%

epvec = 2.^(-8:0.5:-1);
eperr = NaN(numel(epvec), 3);
speed = NaN(size(epvec));

if exist('potstr', 'var') == 0
  potstr = 'log'; 
end

assert(ischar(potstr), 'potstr must be a string');

if exist('N', 'var') == 0
  N = 12e3; 
end

XY = 2 * rand(N, 2) - 1;
V = rand(N, 1);

if exist('maxleaf', 'var') == 0
  maxleaf = 8; 
end

if exist('epk', 'var') == 0
  epk = 1.0e-4; 
end

t0 = tic;
Wref = qtreehat(XY(:, 1), XY(:, 2), V, maxleaf, 0.0, epk, potstr);
t0 = toc(t0);

fprintf(1, 'brute force calc took %f seconds for N = %i random pts\n', t0, N);

for e = 1:numel(epvec)
  t1 = tic;
  Wapx = qtreehat(XY(:, 1), XY(:, 2), V, maxleaf, epvec(e), epk, potstr);
  t1 = toc(t1);
  speed(e) = t1 / t0;

  for c = 1:size(Wref, 2)
    errc = max(abs(Wapx(:, c) - Wref(:, c))) / mean(abs(Wref(:, c))); 
    % fprintf(1, 'rel err:%i = %e @ ep = %f\n', c, errc, epvec(e));
    eperr(e, c) = errc;
  end

  fprintf(1, 'speed ratio = %f @ ep = %f --> err = %e\n', speed(e), epvec(e), eperr(e, 1));
end

if nargout == 0
  % this final (fast) call prints out stats from the MEX program (to be inspected)
  qtreehat(XY(:, 1), XY(:, 2), V, maxleaf, epvec(e), epk, potstr);

  figure; 
  plot(log2(epvec), 1./speed, 'LineWidth', 2, 'Marker', 's');
  line(log2(epvec), ones(size(epvec)), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
  line(log2(epvec), 10 * ones(size(epvec)), 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 2);
  line(log2(epvec), 50 * ones(size(epvec)), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
  grid on;
  xlabel('acc. param. log2(ep)', 'FontSize', 20);
  ylabel('speed-up (compared to brute force)', 'FontSize', 20);
  title(sprintf('N=%i, maxleaf=%i, potstr=%s', N, maxleaf, potstr), 'FontSize', 20);

  figure; 
  plot(log2(epvec), log10(eperr), 'LineWidth', 2, 'Marker', 's');
  line(log2(epvec), -4 * ones(size(epvec)), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
  line(log2(epvec), -3 * ones(size(epvec)), 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 2);
  grid on;
  xlabel('acc. param. log2(ep)', 'FontSize', 20);
  ylabel('log10(error)', 'FontSize', 20);
  hl = legend('value', 'gradx', 'grady');
  set(hl, 'FontSize', 20);
  title(sprintf('N=%i, maxleaf=%i, potstr=%s', N, maxleaf, potstr), 'FontSize', 20);
end

end 
