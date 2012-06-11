
function [dev mean_dev x] = fiber_deviation(est, gt)

for i=1:size(est,2)
    est(:,i) = est(:,i)/norm(est(:,i));
end

for i=1:size(gt,2)
    gt(:,i) = gt(:,i)/norm(gt(:,i));
end

n=size(gt,2);
m=size(est,2);

for i=1:n
    for j=1:m
        dist_mat(i,j)=acos(abs(gt(:,i)'*est(:,j)))*180/pi;
    end
end

[x y]=assignmentoptimal(dist_mat);

mean_dev=y/n;

for k=1:n
 if x(k)~=0
     dev(k,1)=dist_mat(k,x(k));
 end
end

dev=sort(dev);

%        
% 
% deviation=[acos(abs(gt(:,3,1)'*est{1}(:,1)))*180/pi,acos(abs(gt(:,3,1)'*est{1}(:,2)))*180/pi,...
%                      acos(abs(gt(:,3,2)'*est{1}(:,1)))*180/pi,acos(abs(gt(:,3,2)'*est{1}(:,2)))*180/pi];
%                  
% [sorted idx]=sort(deviation);                 
%                  
% if idx(1)==1
%         dev(1)=acos(abs(gt(:,3,1)'*est{1}(:,1)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,2)'*est{1}(:,2)))*180/pi;
%         
% elseif idx(1)==2
%         dev(1)=acos(abs(gt(:,3,1)'*est{1}(:,2)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,2)'*est{1}(:,1)))*180/pi;
%         
% elseif idx(1)==3
%         dev(1)=acos(abs(gt(:,3,2)'*est{1}(:,1)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,1)'*est{1}(:,2)))*180/pi;
%         
% else 
%         dev(1)=acos(abs(gt(:,3,2)'*est{1}(:,2)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,1)'*est{1}(:,1)))*180/pi;
%         
% end
% 
% 
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% end

