close all
clear all
clc

method = 'bas';
est = 'weights_diffus';
weight = 0.2;
snr = 40;

path = sprintf('~/Study/RunningExp/rician_em_weight_est/results/synthetic/%s/%s/', method, est);


d0s = [2 3 4 5]*1.0e-4;
d1s = [1.3:0.2:1.9]*1.0e-3;
angles = [30 60 90];

cc=lines(length(d0s));
h=figure('Position', [100 100 1600 500]);
title(['w_1 = ', weight]);

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
            fileName = sprintf('%s/n=2__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt', path, snr, angles(iAngle), d1, d0, weight);
            dat = dlmread(fileName, '\t');
            meanDev = sum(dat(:,1:2),2)/2;
            maxDev = max(dat(:,1:2),[],2);
            dev = [dev; [mean(meanDev), mean(maxDev)]];
            devstd = [devstd; [std(meanDev), std(maxDev)]];
        end
        figure (h);
        subplot(1,length(angles),iAngle);
        hold all; 
        devHandleBuf(i) = errorbar(d1s', dev(:,1), devstd(:,1), 'color', cc(i,:));
        errorbar(d1s', dev(:,2), devstd(:,2), '--', 'color', cc(i,:));
        hold off;
    end
    
    subplot(1,length(angles), iAngle);
    
    for iLine=1:length(d0s)
        M(iLine,:) = sprintf('d0=%0.1e', d0s(iLine));
    end
    legend(devHandleBuf, M);
    xlabel('Stick diffusivity');
    ylabel('Direction deviation (degrees)');
    originalyLim = ylim;
    disp(originalyLim);
    ylim([originalyLim(1) originalyLim(2)*1.15]);
end

print(h, sprintf('~/Study/HardiToolbox/Summary/figures/synth_weight__snr=%d__w=%0.1f.eps', snr, weight), 'depsc');