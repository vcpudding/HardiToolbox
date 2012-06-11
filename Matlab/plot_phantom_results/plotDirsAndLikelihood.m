clear all
close all
clc

%path = '~/Study/Backup/RunningExp/test_rician_em_diffus/results/phantom/diffus/';
path = '~/Study/Backup/RunningExp/test_rician_em_weight/results/phantom/weights/';

%%%%%plot dirs%%%%%%%%
resFiles = dir([path, 'res*']);
posBuf = [];
dirBuf = [];
hDir = figure;
axis([1 64 1 64]);
axis square;
view(2);
for i=1:length(resFiles)
    
    fileName = resFiles(i).name;
    %disp(fileName);
    pos = sscanf(fileName, 'res_%d_%d_%d.txt');
    
    res = dlmread([path, fileName]);
    nFibers = size(res,2)-1;
    nRep = floor(size(res,1)/5);
    idx = randi(nRep,1,1);
    fibDirs = res(5*(idx-1)+1:5*(idx-1)+3, 1:nFibers);
    posBuf = [posBuf, pos(:, ones(1,nFibers))];
    dirBuf = [dirBuf, fibDirs];
    for j=1:nFibers
        hold on;
        quiver3(pos(1), pos(2), pos(3), fibDirs(1,j), fibDirs(2,j), fibDirs(3,j), 0.3, 'Color', abs(fibDirs(:,j)), 'ShowArrowHead', 'off');
        quiver3(pos(1), pos(2), pos(3), -fibDirs(1,j), -fibDirs(2,j), -fibDirs(3,j), 0.3, 'Color', abs(fibDirs(:,j)), 'ShowArrowHead', 'off');
        hold off;
    end
end


%%%%%%plot likelihood%%%%%%%%%
likeDat = dlmread([path, 'likelihood.txt']');
meanLike = -2000*ones(64);
stdLike = zeros(64);

for i=1:size(likeDat,1)
    dat = likeDat(i,:);
    meanLike(65-dat(2), dat(1)) = dat(4);
    stdLike(65-dat(2), dat(1)) = dat(5);
end

hLike = figure;
subplot(1,2,1); imagesc(meanLike); colorbar; axis square;
subplot(1,2,2); imagesc(stdLike); colorbar; axis square;