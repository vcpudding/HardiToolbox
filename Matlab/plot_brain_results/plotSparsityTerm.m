clear all
close all


x=61;
z=47;
voxType = 'cross';

singleSparseFactor = [1 1e1 1e2 1e3 3e3 5e3 8e3 1e4 1.5e4 2e4 2.5e4 3e4 5e4 8e4 1e5 1.5e5 2e5 2.5e5 3e5 5e5];
crossingSparseFactor = [1 1e1 1e2 1e3 1e4 2e4 4e4 6e4 8e4 1e5 1.5e5 2e5 2.5e5 3e5];
if strcmp(voxType, 'single')
    sparseFactor = singleSparseFactor;
else
    sparseFactor = crossingSparseFactor;
end
weights = [];
sparsity = [];

for s = sparseFactor
    fileName = sprintf('~/Study/RunningExp/rician_em_brain/results/brain/bas/weights_sparse/res_%d_40_%d_3__s=%0.1e.txt',...
        x, z, s);
    res = dlmread(fileName);
    weights = [weights; res(10:13)]; 
    sparsity = [sparsity; sum(res(10:13).^2)];   
end

h = figure;
semilogx(sparseFactor, weights);
hold all;
semilogx(sparseFactor, sparsity, '--');
hold off;
legend('stick 1', 'stick 2', 'stick 3', 'ball', 'sparsity', 'Location', 'NorthWest');
xlabel('Coupling factor');
ylabel('Component weights');

print(h, sprintf('~/Study/HardiToolbox/Summary/figures/sparsity_%d_%d_%s.eps', x, z, voxType), '-depsc');

% 
% for s = crossingSparsity
%     fileName = sprintf('~/Study/RunningExp/rician_em_brain/results/brain/bas/weights_sparse/res_%d_40_%d_3__s=%d.txt',...
%         61, 49, s);
%     res = dlmread(fileName);
%     weights = [weights; res(10:13)];
%     sparsity = [sparsity; sum(res(10:13).^2)];
% end
% 
% figure, semilogx(crossingSparsity, weights);
% hold all;
% semilogx(crossingSparsity, sparsity, '--');
% hold off;
