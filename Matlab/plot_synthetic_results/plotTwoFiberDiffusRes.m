close all
clear all
clc

method = 'modbas';
est = 'weights_diffus';
snr = 40;

path = sprintf('~/Study/RunningExp/rician_em_weight_est/results/synthetic/%s/%s/', method, est);


d0s = [2 3 4 5]*1.0e-4;
d1s = [1.3:0.2:2.1]*1.0e-3;
angles = [30 60 90];
if strcmp(method, 'bas')
    weights = [0.2 0.3 0.4];
else
    weights = [0.3 0.4 0.5];
end

cc=lines(length(d0s));

for iWeight = 1:length(weights);
    weight = weights(iWeight);
for iAngle = 1:length(angles)
    
    h=figure('Position', [100 100 400 400]);
    angle = angles(iAngle);
    devHandleBuf = zeros(length(d0s),1);
    likeBuf = zeros(length(d0s));   
    for i=1:length(d0s)
        d0 = d0s(i);
        dev = [];
        devstd = [];
        like = [];
        likestd = [];
        for d1 = d1s
            fileName = sprintf('%s/n=2__s=%d__a=%d__d1=%0.1e__d0=%0.1e__w=%0.1f.txt', path, snr, angle, d1, d0, weight);
            dat = dlmread(fileName, '\t');
            meanDev = sum(dat(:,1:2),2)/2;
            maxDev = max(dat(:,1:2),[],2);
            dev = [dev; [mean(meanDev), mean(maxDev)]];
            devstd = [devstd; [std(meanDev), std(maxDev)]];
        end
        figure (h);
        %subplot(length(weights),length(angles),(iWeight-1)*length(angles)+iAngle);
        hold all; 
        devHandleBuf(i) = errorbar(d1s', dev(:,1), devstd(:,1), 'color', cc(i,:));
        errorbar(d1s', dev(:,2), devstd(:,2), '--', 'color', cc(i,:));
        hold off;
    end
    
    %subplot(length(weights),length(angles),(iWeight-1)*length(angles)+iAngle);
    if iWeight==1 && iAngle==1
        for iLine=1:length(d0s)
            M(iLine,:) = sprintf('d0=%0.1e', d0s(iLine));
        end
        legend(devHandleBuf, M);
        originalyLim = ylim;
        ylim([originalyLim(1) originalyLim(2)*1.2]);
    end
    xlabel('Stick diffusivity');
    ylabel('Direction deviation (degrees)');
    title(['w_1=', num2str(weight), ', separation angle=', num2str(angle), '\circ']);

    print(h, sprintf('~/Study/HardiToolbox/Summary/figures/synth_%s_%s__snr=%d__w1=%d__angle=%d', method, est, snr, weight*10, angle), '-depsc');
    close(h);
end
end