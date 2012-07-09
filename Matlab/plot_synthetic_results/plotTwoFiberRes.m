close all
clear all
clc

method = 'bas';
est = 'weights_diffus';

path = sprintf('~/Study/RunningExp/rician_em_weight_est/results/synthetic/%s/%s/', method, est);


d0s = [2 3 4 5]*1.0e-4;
d1s = [1.3:0.2:1.9]*1.0e-3;
if strcmp(method, 'bas')
    weights = [0.2 0.3 0.4];
else
    weights = [0.3 0.4 0.5];
end
snr = 40;
% angle = 90;
angles = [30 60 90];


cc=lines(length(d0s));
h=figure('Position', [100 100 1200 500]);

for iAngle = 1:length(angles)
    devHandleBuf = zeros(length(d0s),1);
    likeBuf = zeros(length(d0s));   
for i=1:length(d0s)
    d0 = d0s(i);
    dev = [];
    devstd = [];
    like = [];
    likestd = [];
    for d1 = d1s
        weight = 0.2;
        fileName = sprintf('%s/n=2__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt', path, snr, angles(iAngle), d1, d0, weight);
        dat = dlmread(fileName, '\t');
        meanDev = sum(dat(:,1:2),2)/2;
        maxDev = max(dat(:,1:2),[],2);
        dev = [dev; [mean(meanDev), mean(maxDev)]];
        devstd = [devstd; [std(meanDev), std(maxDev)]];
%         like = [like, mean(dat(:,9))];
%         likestd = [likestd, std(dat(:,9))];
    end
    figure (h);
    subplot(1,length(angles),iAngle);
    hold all; 
    devHandleBuf(i) = errorbar(d1s', dev(:,1), devstd(:,1), 'color', cc(i,:));
    errorbar(d1s', dev(:,2), devstd(:,2), '--', 'color', cc(i,:));
    for iLine=1:length(d0s)
        M(iLine,:) = sprintf('d0=%0.1e', d0s(i));
    end
    legend(devHandleBuf, M);
%     hold off;   
%     subplot(1,2,2);
%     hold on; 
%     likeBuf(i) = errorbar(d1s, like, likestd, 'color', cc(i,:));
    hold off;
end
end

figure(h);
for iSubplot = 1:length(angles)
    subplot(1,length(angles),iSubplot);
    for i=1:length(d0s)
        M(i,:) = sprintf('d0=%0.1e', d0s(i));
    end
    legend(devHandleBuf, M);
%     switch iSubplot
%         case 1
%             legend(devHandleBuf, M);
%             title('Direction deviation');
%             xlabel('Stick diffusivity');
%             ylabel('Direction deviation (degrees)');
%         case 2
%             legend(M);
%             title('Likelihood at convergence');
%             xlabel('Stick diffusivity ');
%             ylabel('Likelihood');
%             %ylim([-1500 -430]);
%     end
end