close all;
clear all;

%load Snew

UnitVectors;
GradientOrientations=g([1:81],:);
NewOrientations=g([1:321],:); % 321
g642=BuildSphere(3);

b=1500;
order=4; delta=500;
sigma=0.009;
v1=[1/sqrt(2) 1/sqrt(2) 0];
orientation1=atan2(v1(2),v1(1));

s=SimulateData(b,GradientOrientations,pi/6,orientation1);

%S=ricernd(Snew,sigma)';

S=DWI_to_SH(s,GradientOrientations,NewOrientations,8,0.006,sigma,false,0.05)'; % regularize and over-sample 81->321

C=constructMatrixOfIntegrals(NewOrientations,order,delta);

c=S*C;

G=constructMatrixOfMonomials(NewOrientations,order);

m=zeros(size(NewOrientations,1),1);

echo on
n=15;

cvx_begin
   variable x1(n)
   dual variables y
   minimize(norm(C*x1-S',2))
   subject to
     y : G * x1 > m;
     %x > 0;
cvx_end

echo off

TensorODF=x1;

%Snew=CreateSignal(x1,g81,order,delta);
% 
% figure,plotTensors(TensorODF,1,[321 1 0]);
W=convert_DT4(TensorODF([15 5 1 12 3 10 11 8 7 14 13 9 4 6 2]));
% W=convertD2W_order6(TensorODF);
[P,V,U]=cp_als(tensor(W),2,'tol',1e-6,'init','nvecs');









