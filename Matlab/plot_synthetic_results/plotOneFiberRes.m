close all
clc

method = 'bas';
est = 'weights_diffus';

path = sprintf('~/Study/RunningExp/rician_em_weight_est/results/synthetic/%s/%s/', method, est);

d0s = [2 3 4 5]*1.0e-4;
d1s = [1.3:0.2:2.1]*1.0e-3;
snr = 20;
bFirst = true;
h=figure('Position', [100 100 1200 500]);
for d0 = d0s
    dev = [];
    devstd = [];
    like = [];
    likestd = [];
    for d1 = d1s
        fileName = sprintf('%s/n=1__s=%d__a=0__d1=%0.1e__d0=%0.1e__w=1.0.txt', path, snr, d1, d0);
        dat = dlmread(fileName, '\t');
        dev = [dev, mean(dat(:,1))];
        devstd = [devstd, std(dat(:,1))];
        like = [like, mean(dat(:,6))];
        likestd = [likestd, std(dat(:,6))];
    end
    figure (h);
    subplot(1,2,1);
    hold all; 
    errorbar(d1s, dev, devstd);  
    hold off;   
    subplot(1,2,2);
    hold all; 
    errorbar(d1s, like, likestd);  
    hold off;
end

figure(h);
for iSubplot = 1:2
    subplot(1,2,iSubplot);
    for i=1:length(d0s)
        M(i,:) = sprintf('d0=%0.1e', d0s(i));
    end
    legend(M);
    switch iSubplot
        case 1
            title('Direction deviation');
            xlabel('Stick diffusivity');
            ylabel('Direction deviation (degrees)');
        case 2
            title('Likelihood at convergence');
            xlabel('Stick diffusivity ');
            ylabel('Likelihood');
    end
end