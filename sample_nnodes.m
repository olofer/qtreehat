function NS = sample_nnodes(N, maxleaf, S, distribstr)
%
% function NS = sample_nnodes(N, maxleaf, S, distribstr)
%
% EXAMPLES:
%   sample_nnodes(1e3, [1,2,4,8,16], 1e3);
%   sample_nnodes(1e3, 8, 1e3, 'uniform');
%

FSZ = 20;

if nargin < 4
  distribstr = 'uniform';
end

if numel(maxleaf) > 1 
  allSamples = cell(numel(maxleaf), 1);
  for l = 1:numel(maxleaf)
    NSl = sample_nnodes(N, maxleaf(l), S, distribstr);
    allSamples{l} = NSl;
    fprintf(1, 'maxleaf = %i\n', maxleaf(l));
  end

  figure;
  lgndStr = cell(numel(maxleaf), 1);
  hold on;
  for l = 1:numel(maxleaf)
    NSl = allSamples{l};
    plot(NSl(:, 1), NSl(:, 2), 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 8);    
    lgndStr{l} = sprintf('maxleaf = %i', maxleaf(l));
  end
  grid on;
  xlabel('number of points', 'FontSize', FSZ);
  ylabel('number of nodes', 'FontSize', FSZ);
  hl = legend(lgndStr);
  set(hl, 'FontSize', FSZ);
  title(sprintf('%i (%s) samples each', S, distribstr), 'FontSize', FSZ);

  NS = allSamples;
  return;
end

assert(numel(maxleaf) == 1);

NS = NaN(S, 2); % columns = {N, nno}

for s = 1:S 
  Ns = max([round(rand * N), maxleaf]);
  if strcmp(distribstr, 'uniform')
    XY = 2 * rand(Ns, 2) - 1;
  elseif strcmp(distribstr, 'normal')
    XY = randn(Ns, 2);
  elseif strcmp(distribstr, 'twogroups')
    c1 = floor(Ns / 2);
    c2 = Ns - c1;
    if c1 >= 1 && c2 >= 1
      XY = 0.1 * randn(c1, 2) + repmat(rand(1, 2), [c1, 1]);
      XY = [XY; 0.1 * randn(c2, 2) + repmat(rand(1, 2), [c2, 1])]; 
    else
      XY = randn(Ns, 2);
    end
  end
  assert(size(XY, 1) == Ns);
  assert(size(XY, 2) == 2);
  V = rand(Ns, 1);

  %[~, nno] = qtreehat(XY(:, 1), XY(:, 2), V, maxleaf, 0.50, 0.0, 'log');
  
  [dummy, nno] = qtreehat(XY(:, 1), XY(:, 2), V, maxleaf, 1e-6, 0.0, 'log', [0.0], [0.0]);
  assert(size(dummy, 1) == 1 && size(dummy, 2) == 3);

  NS(s, :) = [Ns, nno];
end

if nargout == 1
  return;
end

NB = NaN(N, 1);
for b = 1:N 
  NB(b) = bound_nodes(b, maxleaf);
end

figure;
plot(NS(:, 1), NS(:, 2), 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 8);
hold on;
%plot(1:N, ((1:N) * (1 + maxleaf)) / maxleaf, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
plot(1:N, NB, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 2);
grid on;
xlabel('number of points', 'FontSize', FSZ);
ylabel('number of nodes', 'FontSize', FSZ);
title(sprintf('%i (%s) samples, using maxleaf = %i', S, distribstr, maxleaf), 'FontSize', FSZ);

end

function NN = bound_nodes(N, L)
  % this is the same expression as used in 'qtreehat.c'; more tight bound would be preferred
  NN = min([16 / 3, L + 1]) * N/L + 128;
end
