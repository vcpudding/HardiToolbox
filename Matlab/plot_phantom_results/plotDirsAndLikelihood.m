clear all
close all
clc

method = 'bas';
est = 'diffus';

%path = '~/Study/Backup/RunningExp/test_rician_em_diffus/results/phantom/diffus/';
%path = '~/Study/RunningExp/rician_em_phantom/results/phantom/bas/weights/';
path = sprintf('~/Study/RunningExp/rician_em_phantom/results/phantom/%s/%s/', method, est);

%%%%%plot dirs%%%%%%%%
resFiles = dir([path, 'res*']);
hDir = figure('Position', [200 200 400 400]);
axis([5 60 5 55]);
axis square;
view(2);
meanLike = ones(64);
stdLike = zeros(64);
for i=1:length(resFiles)
    
    fileName = resFiles(i).name;
    %disp(fileName);
    pos = sscanf(fileName, 'res_%d_%d_%d_%d.txt');
    
    res = dlmread([path, fileName]);
    if size(res,2)<7
        disp(fileName);
        continue;
    end
    if strcmp(method, 'bas')
        nFibers = (size(res,2)-3)/5;
    else
        nFibers = (size(res,2)-2)/5;
    end
    nRep = size(res,1);
    idx = randi(nRep,1,1);
    fibDirs = res(idx, 1:3*nFibers);
    fibDirs = reshape(fibDirs, [3, nFibers]);
    for j=1:nFibers
        hold on;
        quiver3(pos(1), pos(2), pos(3), fibDirs(1,j), fibDirs(2,j), fibDirs(3,j), 0.3, 'Color', abs(fibDirs(:,j)), 'ShowArrowHead', 'off');
        quiver3(pos(1), pos(2), pos(3), -fibDirs(1,j), -fibDirs(2,j), -fibDirs(3,j), 0.3, 'Color', abs(fibDirs(:,j)), 'ShowArrowHead', 'off');
        hold off;
    end
    
    likeDat = res(:, size(res,2));
    meanLike(65-pos(2), pos(1)) = mean(likeDat);
    stdLike(65-pos(2), pos(1)) = std(likeDat);
end

meanLike(meanLike>0) = min(meanLike(:))-100;

%%%%%%plot likelihood%%%%%%%%%

hLike = figure('Position', [100 100 1200 400]);
subplot(1,2,1); imagesc(meanLike); colorbar; axis square;
subplot(1,2,2); imagesc(stdLike); colorbar; axis square;

print(hDir, sprintf('~/Study/HardiToolbox/Summary/figures/phantom_%s_%s_dir.eps', method, est), '-depsc');
print(hLike, sprintf('~/Study/HardiToolbox/Summary/figures/phantom_%s_%s_like.eps', method, est), '-depsc');