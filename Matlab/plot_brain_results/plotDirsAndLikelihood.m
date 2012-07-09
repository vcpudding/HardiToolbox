clear all
close all
clc

method = 'bas';
est = 'weights';

path = sprintf('~/Study/RunningExp/rician_em_brain/results/brain/%s/%s/', method, est);
load('~/Study/Data/Brain_data/FA_SLICE40_CORONAL_B2000.mat');

%%%%%plot dirs%%%%%%%%
resFiles = dir([path, 'res*_3.txt']);
hDir = figure('Position', [200 200 400 400]);
axis([0 22 0 22]);
axis square;
view(2);
meanLike = ones(21);
stdLike = zeros(21);
weightThres = 0.10;
for i=1:length(resFiles)
    
    fileName = resFiles(i).name;
    pos = sscanf(fileName, 'res_%d_%d_%d_%d__s=1e5.txt');
    
%     if pos(1)==65 && pos(3)==59
%         stop = 1;
%     end
    if pos(1)==52 && pos(3)==44
        stop=1;
    end
    
    pos(1) = pos(1)-49;
    pos(2) = 0;
    pos(3) = pos(3)-39;
    
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
    f = fa(pos(1), pos(3));
    weights = res(idx, 3*nFibers+1:3*nFibers+nFibers);
    fibDirs = fibDirs(:, weights>weightThres);
    weights = weights(weights>weightThres);
    %weights = weights/(sum(weights)+1e-6);
    nFibers = length(weights);
    scale = weights;
    if max(scale)<1e-3
        continue;
    end
    for j=1:nFibers
        hold on;
        quiver3(pos(1), pos(3), 0, fibDirs(1,j)*scale(j), fibDirs(3,j)*scale(j), fibDirs(2,j)*scale(j), 1, 'Color', abs(fibDirs(:,j)), 'ShowArrowHead', 'off');
        quiver3(pos(1), pos(3), 0, -fibDirs(1,j)*scale(j), -fibDirs(3,j)*scale(j), -fibDirs(2,j)*scale(j), 1, 'Color', abs(fibDirs(:,j)), 'ShowArrowHead', 'off');
        hold off;
    end
    
    likeDat = res(:, size(res,2));
    meanLike(22-pos(3), pos(1)) = mean(likeDat);
    stdLike(22-pos(3), pos(1)) = std(likeDat);
end

meanLike(meanLike>0) = min(meanLike(:))-100;

%%%%%%plot likelihood%%%%%%%%%

% hLike = figure('Position', [100 100 1200 400]);
% subplot(1,2,1); imagesc(meanLike); colorbar; axis square;
% subplot(1,2,2); imagesc(stdLike); colorbar; axis square;
% 
% print(hDir, sprintf('~/Study/HardiToolbox/Summary/figures/brain_%s_%s_dir.eps', method, est), '-depsc');
% print(hLike, sprintf('~/Study/HardiToolbox/Summary/figures/brain_%s_%s_like.eps', method, est), '-depsc');