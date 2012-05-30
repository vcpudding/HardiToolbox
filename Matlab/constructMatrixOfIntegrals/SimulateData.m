function S=SimulateData(b,GradientDirections,angle,orientation1)

l=[1.7e-3, 3e-4, 3e-4];

R=[cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 01];
%R1=[1/sqrt(2) 1/sqrt(2) 0;-1/sqrt(2) 1/sqrt(2) 0;0 0 1]; % fiber orientation 1
R1=[cos(orientation1) sin(orientation1) 0;-sin(orientation1) cos(orientation1) 0;0 0 1];
R2=R1*R;
R3=R2*R;

L=diag(l);

D1=R1'*L*R1;
D2=R2'*L*R2;
D3=R3'*L*R3;

p1=0.5;p2=0.5; p3=0; % for two fibers
%p1=0.4; p2=0.3; p3=0.3; % for three fibers

for i=1:size(GradientDirections,1);
    S(i)=p1*exp(-b*GradientDirections(i,:)*D1*GradientDirections(i,:)')+p2*exp(-b*GradientDirections(i,:)*D2*GradientDirections(i,:)')+p3*exp(-b*GradientDirections(i,:)*D3*GradientDirections(i,:)');
end



    

