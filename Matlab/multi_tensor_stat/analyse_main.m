
clear all
close all

angles=[30 45 60 90];

for k=1:length(angles)

% f1='n=2__b=3000__g=64__a=';
% f2=num2str(angles(k));
% f3='__s=40__w=0.5_d.txt';

%filename=strcat(f1,f2,f3);
filename = sprintf('../../results/n=2__est=2__b=3000__g=64__a=%d__s=%d__w=0.5__init=1_d.txt', angles(k), 20000);

Directions=load(filename);

%Directions=n_2__b_3000__g_64__a_60__s_20__w_0_5_d;

no_exp=size(Directions,1)/3;

for i=1:no_exp
    Est_dir(:,1)=Directions((i-1)*3+1:3*i,3);
    Est_dir(:,2)=Directions((i-1)*3+1:3*i,4);
    True_dir(:,1)=Directions((i-1)*3+1:3*i,1);
    True_dir(:,2)=Directions((i-1)*3+1:3*i,2);
    
   [dev mean_dev x] = fiber_deviation(Est_dir,True_dir);
   
   Dev(i,:)=dev;
   
   Mean_Dev(i,k)=mean_dev;
   
end

end

figure,errorbar(angles, mean(Mean_Dev),std(Mean_Dev),'.');

   
   
    

